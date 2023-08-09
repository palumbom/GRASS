function precompute_quantities_gpu!(disk::DiskParams{T}, gpu_allocs::GPUAllocs{T}) where T<:AF
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
        θc_subgrid[idx, 1:num_edges-1] .= GRASS.get_grid_centers(edges)'
        θe_subgrid[idx, 1:num_edges] .= collect(edges)'
    end

    # copy data to GPU
    @cusync begin
        # get observer vectoir and rotation matrix
        O⃗ = CuArray{precision}(disk.O⃗)
        R_θ = CuArray{precision}(disk.R_θ)

        # move grids to GPU
        ϕe_subgrid_gpu = CuArray{precision}(ϕe_subgrid)
        ϕc_subgrid_gpu = CuArray{precision}(ϕc_subgrid)
        θe_subgrid_gpu = CuArray{precision}(θe_subgrid)
        θc_subgrid_gpu = CuArray{precision}(θc_subgrid)

        # allocate memory for computations on subgrid
        vec1 = CUDA.zeros(precision, size(θc_subgrid)..., 3)
        vec2 = CUDA.zeros(precision, size(θc_subgrid)..., 3)
        vec3 = CUDA.zeros(precision, size(θc_subgrid)..., 3)

        # allocat larger arrays
        μs = CUDA.zeros(precision, size(θc_subgrid))
        wts = CUDA.zeros(precision, size(θc_subgrid))
        z_rot = CUDA.zeros(precision, size(θc_subgrid))
    end

    # alias from GPU allocs
    ax_codes = gpu_allocs.ax_codes

    # compute quantities at subgrid points
    threads1 = (16,16)
    blocks1 = cld(prod(size(θe_subgrid)), prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 precompute_quantities_gpu!(vec1, vec2, vec3, μs, wts, z_rot,
                                                                                       ϕe_subgrid_gpu, θe_subgrid_gpu,
                                                                                       ϕc_subgrid_gpu, θc_subgrid_gpu,
                                                                                       R_θ, O⃗, ρs, A, B, C, u1, u2)

    # average quantities over subgrid
    threads1 = 256
    blocks1 = cld(disk.Nsubgrid^2, prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 average_subgrid_gpu!(μs, wts, z_rot, vec1, ax_codes, disk.N)

    # copy over to gpu allocs
    idx = CartesianIndices((1:disk.N, 1:disk.N))
    @cusync begin
        CUDA.copyto!(gpu_allocs.μs, idx, μs, idx)
        CUDA.copyto!(gpu_allocs.wts, idx, wts, idx)
        CUDA.copyto!(gpu_allocs.z_rot, idx, z_rot, idx)
    end

    # instruct CUDA to free up unneeded memory
    @cusync begin
        CUDA.unsafe_free!(O⃗)
        CUDA.unsafe_free!(R_θ)
        CUDA.unsafe_free!(ϕe_subgrid_gpu)
        CUDA.unsafe_free!(ϕc_subgrid_gpu)
        CUDA.unsafe_free!(θe_subgrid_gpu)
        CUDA.unsafe_free!(θc_subgrid_gpu)
        CUDA.unsafe_free!(vec1)
        CUDA.unsafe_free!(vec2)
        CUDA.unsafe_free!(vec3)
        CUDA.unsafe_free!(μs)
        CUDA.unsafe_free!(wts)
        CUDA.unsafe_free!(z_rot)
    end

    return nothing
end

function precompute_quantities_gpu!(vec1, vec2, vec3, μs, wts, z_rot, ϕe, θe, ϕc, θc, R_θ, O⃗, ρs, A, B, C, u1, u2)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # loop over subtiled grid
    for i in idx:sdx:CUDA.length(ϕc)
        for j in idy:sdy:CUDA.size(θc,1)
            # move on if we are done with the tiles in latitude slice
            if j > 1 && θc[i,j] == 0.0
                continue
            end

            # take views of pre-allocated memory
            xyz = CUDA.view(vec1, i, j, :)
            abc = CUDA.view(vec2, i, j, :)
            def = CUDA.view(vec3, i, j, :)

            # get cartesian coords
            sphere_to_cart_gpu!(xyz, ρs, ϕc[i], θc[i,j])

            # get vector from spherical circle center to surface patch
            @inbounds abc[1] = CUDA.copy(xyz[1])
            @inbounds abc[2] = CUDA.copy(xyz[2])
            @inbounds abc[3] = 0.0

            # take cross product to get vector in direction of rotation
            @inbounds def[1] = CUDA.copy(abc[2] * ρs)
            @inbounds def[2] = CUDA.copy(- abc[1] * ρs)
            @inbounds def[3] = 0.0

            # make it a unit vector
            def_norm = CUDA.sqrt(def[1]^2.0 + def[2]^2.0)
            @inbounds def[1] /= def_norm
            @inbounds def[2] /= def_norm

            # set magnitude by differential rotation
            v0 = 0.000168710673 # Rsol/day/speed of light
            rp = (v0 / rotation_period_gpu(ϕc[i], A, B, C))
            @inbounds def[1] *= rp
            @inbounds def[2] *= rp

            # rotate xyz by inclination and calculate mu
            rotate_vector_gpu!(xyz, R_θ)
            @inbounds μs[i,j] = calc_mu_gpu(xyz, O⃗)
            if μs[i,j] <= 0.0
                continue
            end

            # rotate the velocity vectors by inclination
            rotate_vector_gpu!(def, R_θ)

            # get vector pointing from observer to surface patch
            abc[1] = CUDA.copy(xyz[1] - O⃗[1])
            abc[2] = CUDA.copy(xyz[2] - O⃗[2])
            abc[3] = CUDA.copy(xyz[3] - O⃗[3])

            # get angle between them
            n1 = CUDA.sqrt(abc[1]^2.0 + abc[2]^2.0 + abc[3]^2.0)
            n2 = CUDA.sqrt(def[1]^2.0 + def[2]^2.0 + def[3]^2.0)
            angle = (abc[1] * def[1] + abc[2] * def[2] + abc[3] * def[3])
            angle /= (n1 * n2)

            # project velocity onto line of sight
            @inbounds z_rot[i,j] = n2 * angle

            # calculate the limb darkening
            ld = quad_limb_darkening(μs[i,j], u1, u2)

            # get projected area element
            dϕ = ϕe[i+1] - ϕe[i]
            dθ = θe[i,j+1] - θe[i,j]
            dA = calc_dA(ρs, ϕc[i], dϕ, dθ)
            @inbounds abc[1] = CUDA.copy(xyz[1] - O⃗[1])
            @inbounds abc[2] = CUDA.copy(xyz[2] - O⃗[2])
            @inbounds abc[3] = CUDA.copy(xyz[3] - O⃗[3])
            dA *= CUDA.abs(abc[1] * xyz[1] + abc[2] * xyz[2] + abc[3] * xyz[3])

            # get weights as product of limb darkening and projected dA
            @inbounds wts[i,j] = ld * dA
        end
    end
    return nothing
end

function average_subgrid_gpu!(μs, wts, z_rot, xyz, ax_codes, N)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    # get number of elements along tile side
    k = Int(CUDA.size(μs,1) / N)

    # total number of elements output array
    num_tiles = N^2

    for itr in idx:sdx:num_tiles
        # get index for output array
        row = (itr - 1) ÷ N
        col = (itr - 1) % N

        # get indices for input array
        i = row * k + 1
        j = col * k + 1

        # set up sum holders for scalars
        μ_sum = CUDA.zero(CUDA.eltype(μs))
        w_sum = CUDA.zero(CUDA.eltype(wts))
        v_sum = CUDA.zero(CUDA.eltype(z_rot))

        # set up sum holders for vector components
        x_sum = CUDA.zero(CUDA.eltype(xyz))
        y_sum = CUDA.zero(CUDA.eltype(xyz))
        z_sum = CUDA.zero(CUDA.eltype(xyz))

        # initiate counter
        count = 0

        for ti in i:i+k-1, tj in j:j+k-1
            if μs[ti, tj] > 0
                # sum on scalar quantities
                μ_sum += μs[ti, tj]
                w_sum += wts[ti, tj]
                v_sum += z_rot[ti, tj]

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
            @inbounds μs[row + 1, col + 1] = μ_sum / count
            @inbounds wts[row + 1, col + 1] = w_sum / count
            @inbounds z_rot[row + 1, col + 1] = v_sum / count

            # set scalar quantity elements as average
            @inbounds xyz[row + 1, col + 1, 1] = x_sum / count
            @inbounds xyz[row + 1, col + 1, 2] = y_sum / count
            @inbounds xyz[row + 1, col + 1, 3] = z_sum / count
        else
            # zero out elements if count is 0 (avoid div by 0)
            @inbounds μs[row + 1, col + 1] = 0.0
            @inbounds wts[row + 1, col + 1] = 0.0
            @inbounds z_rot[row + 1, col + 1] = 0.0

            @inbounds xyz[row + 1, col + 1, 1] = x_sum / count
            @inbounds xyz[row + 1, col + 1, 2] = y_sum / count
            @inbounds xyz[row + 1, col + 1, 3] = z_sum / count
        end

        # while we are here, compute ax codes
        x = xyz[row + 1, col + 1, 1]
        z = xyz[row + 1, col + 1, 3]
        @inbounds ax_codes[row + 1, col + 1] = find_nearest_ax_gpu(x, z)
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
