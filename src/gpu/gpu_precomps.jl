function precompute_quantities_gpu!(disk::DiskParams{T1}, gpu_allocs::GPUAllocs{T2}) where {T1<:AF, T2<:AF}
    # get precision from GPU allocs
    precision = eltype(gpu_allocs.λs)

    # convert scalars from disk params to desired precision
    ρs = convert(precision, disk.ρs)
    A = convert(precision, disk.A)
    B = convert(precision, disk.B)
    C = convert(precision, disk.C)
    u1 = convert(precision, disk.u1)
    u2 = convert(precision, disk.u2)

    # get subtile edges for subgrid'd latitude
    N_edges_subgrid = disk.Nsubgrid * disk.N + 1
    ϕe_subgrid = range(first(disk.ϕe), last(disk.ϕe), length=N_edges_subgrid)
    ϕc_subgrid = get_grid_centers(ϕe_subgrid)

    # allocate for longitude grid
    θc_subgrid = zeros(length(ϕc_subgrid), maximum(disk.Nθ) * disk.Nsubgrid)
    θe_subgrid = zeros(size(θc_subgrid) .+ 1 ...)

    # loop over latitutde slices in original course grid
    for i in eachindex(disk.Nθ)
        # get idx for longitude subgrid slices in courser lat slice
        idx  = (i - 1) * disk.Nsubgrid+1:(i * disk.Nsubgrid)

        # get number of edges needed and make grids
        num_edges = disk.Nθ[i] * disk.Nsubgrid + 1
        edges = range(deg2rad(0.0), deg2rad(360.0), length=num_edges)
        θc_subgrid[idx, 1:num_edges-1] .= get_grid_centers(edges)'
        θe_subgrid[idx, 1:num_edges] .= collect(edges)'
    end

    # get number of tiles in each row of subgrid
    Nθ_sub = map(x -> findfirst(x[2:end] .== 0.0), eachrow(θc_subgrid))
    Nθ_sub[isnothing.(Nθ_sub)] .= size(θc_subgrid, 2)
    Nθ_sub = Array{Int}(Nθ_sub)

    # copy data to GPU
    @cusync begin
        # get observer vectoir and rotation matrix
        O⃗ = CuArray{precision}(disk.O⃗)
        Nθ = CuArray{Int32}(Nθ_sub)
        R_θ = CuArray{precision}(disk.R_θ)

        # move grids to GPU
        ϕe_subgrid_gpu = CuArray{precision}(ϕe_subgrid)
        ϕc_subgrid_gpu = CuArray{precision}(ϕc_subgrid)
        θe_subgrid_gpu = CuArray{precision}(θe_subgrid)
        θc_subgrid_gpu = CuArray{precision}(θc_subgrid)

        # allocate memory for computations on subgrid
        xyz = CUDA.zeros(precision, size(θc_subgrid)..., 3)

        # allocate larger arrays
        μs = CUDA.zeros(precision, size(θc_subgrid))
        ld = CUDA.zeros(precision, size(θc_subgrid))
        dA = CUDA.zeros(precision, size(θc_subgrid))
        z_rot = CUDA.zeros(precision, size(θc_subgrid))
    end

    # alias from GPU allocs
    μs_out = gpu_allocs.μs
    wts_out = gpu_allocs.wts
    z_rot_out = gpu_allocs.z_rot
    ax_codes = gpu_allocs.ax_codes

    # compute quantities at subgrid points
    threads1 = (16,16)
    blocks1 = cld(disk.N^2, prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 precompute_quantities_gpu!(xyz, μs, ld, dA, z_rot,
                                                                                       ϕe_subgrid_gpu, θe_subgrid_gpu,
                                                                                       ϕc_subgrid_gpu, θc_subgrid_gpu,
                                                                                       Nθ, R_θ, O⃗, ρs, A, B, C, u1, u2)

    # average quantities over subgrid
    threads1 = 256
    blocks1 = cld(disk.Nsubgrid^2, prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 average_subgrid_gpu!(μs_out, μs, wts_out, ld, dA,
                                                                                 z_rot_out, z_rot, xyz,
                                                                                 ax_codes, disk.N)

    # instruct CUDA to free up unneeded memory
    @cusync begin
        CUDA.unsafe_free!(O⃗)
        CUDA.unsafe_free!(Nθ)
        CUDA.unsafe_free!(R_θ)
        CUDA.unsafe_free!(ϕe_subgrid_gpu)
        CUDA.unsafe_free!(ϕc_subgrid_gpu)
        CUDA.unsafe_free!(θe_subgrid_gpu)
        CUDA.unsafe_free!(θc_subgrid_gpu)
        CUDA.unsafe_free!(xyz)
        CUDA.unsafe_free!(μs)
        CUDA.unsafe_free!(ld)
        CUDA.unsafe_free!(dA)
        CUDA.unsafe_free!(z_rot)
    end

    # instruct the garbage collect to clean up GPU memory
    GC.gc(true)

    CUDA.synchronize()
    return nothing
end

function precompute_quantities_gpu!(all_xyz, μs, ld, dA, z_rot, ϕe, θe, ϕc, θc, Nθ, R_θ, O⃗, ρs, A, B, C, u1, u2)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # loop over subtiled grid
    for i in idx:sdx:CUDA.length(ϕc)
        for j in idy:sdy:Nθ[i]
            # take views of pre-allocated memory
            xyz = CUDA.view(all_xyz, i, j, :)

            # get cartesian coords
            x, y, z = sphere_to_cart_gpu(ρs, ϕc[i], θc[i,j])

            # get vector from spherical circle center to surface patch
            a = x
            b = y
            c = 0.0

            # take cross product to get vector in direction of rotation
            d = b * ρs
            e = - a * ρs
            f = 0.0

            # make it a unit vector
            def_norm = CUDA.sqrt(d^2.0 + e^2.0)
            d /= def_norm
            e /= def_norm

            # set magnitude by differential rotation
            v0 = Float64(0.000168710673) # Rsol/day/speed of light
            rp = (v0 / rotation_period_gpu(ϕc[i], A, B, C))
            d *= rp
            e *= rp

            # rotate xyz by inclination, store vector for later
            x, y, z = rotate_vector_gpu(x, y, z, R_θ)
            @inbounds xyz[1] = x
            @inbounds xyz[2] = y
            @inbounds xyz[3] = z

            # rotate xyz by inclination and calculate mu
            @inbounds μs[i,j] = calc_mu_gpu(x, y, z, O⃗)
            if μs[i,j] <= 0.0
                continue
            end

            # rotate the velocity vectors by inclination
            d, e, f = rotate_vector_gpu(d, e, f, R_θ)

            # get vector pointing from observer to surface patch
            a = x - O⃗[1]
            b = y - O⃗[2]
            c = z - O⃗[3]

            # get angle between them
            n1 = CUDA.sqrt(a^2.0 + b^2.0 + c^2.0)
            n2 = CUDA.sqrt(d^2.0 + e^2.0 + f^2.0)
            angle = (a * d + b * e + c * f)
            angle /= (n1 * n2)

            # project velocity onto line of sight
            @inbounds z_rot[i,j] = n2 * angle

            # calculate the limb darkening
            ld[i,j] = quad_limb_darkening(μs[i,j], u1, u2)

            # get area element
            dϕ = ϕe[i+1] - ϕe[i]
            dθ = θe[i,j+1] - θe[i,j]
            @inbounds dA[i,j] = calc_dA(ρs, ϕc[i], dϕ, dθ)

            # project onto line of sight
            @inbounds dA[i,j] *= CUDA.abs(a * x + b * y + c * z)
        end
    end
    return nothing
end

function average_subgrid_gpu!(μs_out, μs, wts_out, ld, dA, z_rot_out, z_rot, xyz, ax_codes, N)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    # get number of elements along tile side
    k = Int(CUDA.size(μs,1) / N)

    # total number of elements output array
    num_tiles = N^2

    for t in idx:sdx:num_tiles
        # get index for output array
        row = (t - 1) ÷ N
        col = (t - 1) % N

        # get indices for input array
        i = row * k + 1
        j = col * k + 1

        # set up sum holders for scalars
        μ_sum = CUDA.zero(CUDA.eltype(μs))
        v_sum = CUDA.zero(CUDA.eltype(z_rot))
        ld_sum = CUDA.zero(CUDA.eltype(ld))
        dA_sum = CUDA.zero(CUDA.eltype(dA))

        # set up sum holders for vector components
        x_sum = CUDA.zero(CUDA.eltype(xyz))
        y_sum = CUDA.zero(CUDA.eltype(xyz))
        z_sum = CUDA.zero(CUDA.eltype(xyz))

        # initiate counter
        count = 0

        for ti in i:i+k-1, tj in j:j+k-1
            if μs[ti, tj] > 0.0
                # sum on scalar quantities
                μ_sum += μs[ti, tj]
                v_sum += z_rot[ti, tj]
                ld_sum += ld[ti, tj]
                dA_sum += dA[ti, tj]

                # sum on vector components
                x_sum += xyz[ti, tj, 1]
                y_sum += xyz[ti, tj, 2]
                z_sum += xyz[ti, tj, 3]

                # iterate counter
                count += 1
            end
        end

        if count > 0
            # set scalar quantity elements as average
            @inbounds μs_out[row + 1, col + 1] = μ_sum / count
            @inbounds z_rot_out[row + 1, col + 1] = v_sum / count
            @inbounds wts_out[row + 1, col + 1] = (ld_sum / count) * dA_sum

            # set scalar quantity elements as average
            @inbounds xx = x_sum / count
            @inbounds yy = y_sum / count
            @inbounds zz = z_sum / count

            # get axis code
            @inbounds ax_codes[row + 1, col + 1] = find_nearest_ax_gpu(xx, zz)
        else
            # zero out elements if count is 0 (avoid div by 0)
            @inbounds μs_out[row + 1, col + 1] = 0.0
            @inbounds wts_out[row + 1, col + 1] = 0.0
            @inbounds z_rot_out[row + 1, col + 1] = 0.0

            # set axis code to 0
            @inbounds ax_codes[row + 1, col + 1] = 0
        end
    end
    return nothing
end

function get_keys_and_cbs_gpu!(gpu_allocs::GPUAllocs{T}, soldata::GPUSolarData{T}) where T<:AF
    # parse out gpu allocs
    μs = gpu_allocs.μs
    z_cbs = gpu_allocs.z_cbs
    dat_idx = gpu_allocs.dat_idx
    ax_codes = gpu_allocs.ax_codes

    # parse out soldata
    cbsall = soldata.cbs
    disc_mu = soldata.mu
    disc_ax = soldata.ax

    threads1 = (16, 16)
    blocks1 = cld(prod(size(μs)), prod(threads1))

    @cusync @captured @cuda threads=threads1 blocks=blocks1 get_keys_and_cbs_gpu!(dat_idx, z_cbs, μs, ax_codes,
                                                                                  cbsall, disc_mu, disc_ax)
    CUDA.synchronize()
    return nothing
end

function get_keys_and_cbs_gpu!(dat_idx, z_cbs, μs, ax_codes, cbsall, disc_mu, disc_ax)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    for i in idx:sdx:CUDA.size(μs,1)
        for j in idy:sdy:CUDA.size(μs,2)
            # move to next iter if off disk
            if μs[i,j] <= 0.0
                continue
            end

            # find the data index for location on disk
            idx = find_data_index_gpu(μs[i,j], ax_codes[i,j], disc_mu, disc_ax)
            dat_idx[i,j] = idx
            z_cbs[i,j] = cbsall[idx]
        end
    end
    return nothing
end
