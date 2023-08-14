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

    # get size of sub-tiled grid
    N = disk.N
    Nsubgrid = disk.Nsubgrid
    Nϕ_sub = N * Nsubgrid
    Nθ_sub = maximum(disk.Nθ) * Nsubgrid

    # copy data to GPU
    @cusync begin
        # get observer vectoir and rotation matrix
        O⃗ = CuArray{precision}(disk.O⃗)
        Nθ = CuArray{Int32}(disk.Nθ)
        R_x = CuArray{precision}(disk.R_x)

        # allocate memory for computations on subgrid
        μs = CUDA.zeros(precision, Nϕ_sub, Nθ_sub)
        ld = CUDA.zeros(precision, Nϕ_sub, Nθ_sub)
        dA = CUDA.zeros(precision, Nϕ_sub, Nθ_sub)
        xz = CUDA.zeros(precision, Nϕ_sub, Nθ_sub, 2)
        z_rot = CUDA.zeros(precision, Nϕ_sub, Nθ_sub)
    end

    # compute quantities at subgrid points
    threads1 = (16,16)
    blocks1 = cld(Nsubgrid^2, prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 precompute_quantities_gpu!(xz, μs, ld, dA, z_rot, N, Nsubgrid,
                                                                                       Nθ, R_x, O⃗, ρs, A, B, C, u1, u2)

    # alias from GPU allocs
    μs_out = gpu_allocs.μs
    wts_out = gpu_allocs.wts
    z_rot_out = gpu_allocs.z_rot
    ax_codes = gpu_allocs.ax_codes

    # average quantities over subgrid
    threads1 = 256
    blocks1 = cld(Nsubgrid^2, prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 average_subgrid_gpu!(μs_out, μs, wts_out, ld, dA,
                                                                                 z_rot_out, z_rot, xz,
                                                                                 ax_codes, disk.N, maximum(disk.Nθ))

    # instruct CUDA to free up unneeded memory
    @cusync begin
        CUDA.unsafe_free!(O⃗)
        CUDA.unsafe_free!(Nθ)
        CUDA.unsafe_free!(R_x)
        CUDA.unsafe_free!(μs)
        CUDA.unsafe_free!(ld)
        CUDA.unsafe_free!(dA)
        CUDA.unsafe_free!(xz)
        CUDA.unsafe_free!(z_rot)
    end

    # instruct the garbage collect to clean up GPU memory
    GC.gc(true)

    CUDA.synchronize()
    return nothing
end

function precompute_quantities_gpu!(all_xz, μs, ld, dA, z_rot, N, Nsubgrid, Nθ, R_x, O⃗, ρs, A, B, C, u1, u2)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # get latitude subtile step size
    N_ϕ_edges = N * Nsubgrid
    dϕ = (CUDA.deg2rad(90.0) - CUDA.deg2rad(-90.0)) / (N_ϕ_edges)

    # loop over subtiled latitudes
    for i in idx:sdx:(N_ϕ_edges)
        # get coordinates of latitude subtile center
        ϕc = CUDA.deg2rad(-90.0) + (dϕ/2.0) + (i - 1) * dϕ

        # get number of longitude tiles in course latitude slice
        k = CUDA.div(i - 1, Nsubgrid) + 1
        N_θ_edges = Nθ[k] * Nsubgrid

        # get longitude subtile step size
        dθ = (CUDA.deg2rad(360.0) - CUDA.deg2rad(0.0)) / (N_θ_edges)

        # loop over subtiled latitudes
        for j in idy:sdy:N_θ_edges
            # get longitude
            θc = CUDA.deg2rad(0.0) + (dθ/2.0) + (j - 1) * dθ

            # take views of pre-allocated memory
            xz = CUDA.view(all_xz, i, j, :)

            # get cartesian coords
            x, y, z = sphere_to_cart_gpu(ρs, ϕc, θc)

            # get vector from spherical circle center to surface patch
            a = x
            b = y
            c = CUDA.zero(CUDA.eltype(μs))

            # take cross product to get vector in direction of rotation
            d = b * ρs
            e = - a * ρs
            f = CUDA.zero(CUDA.eltype(μs))

            # make it a unit vector
            def_norm = CUDA.sqrt(d^2.0 + e^2.0)
            d /= def_norm
            e /= def_norm

            # set magnitude by differential rotation
            v0 = 0.000168710673 # Rsol/day/speed of light
            rp = (v0 / rotation_period_gpu(ϕc, A, B, C))
            d *= rp
            e *= rp

            # rotate xyz by inclination
            x, y, z = rotate_vector_gpu(x, y, z, R_x)

            # rotate xyz by inclination and calculate mu
            @inbounds μs[i,j] = calc_mu_gpu(x, y, z, O⃗)
            if μs[i,j] <= 0.0
                continue
            end

            # rotate the velocity vectors by inclination
            d, e, f = rotate_vector_gpu(d, e, f, R_x)

            # get vector pointing from observer to surface patch
            a = x - O⃗[1]
            b = y - O⃗[2]
            c = z - O⃗[3]

            # get angle between them
            n1 = CUDA.sqrt(a^2.0 + b^2.0 + c^2.0)
            n2 = CUDA.sqrt(d^2.0 + e^2.0 + f^2.0)
            angle = (a * d + b * e + c * f) / (n1 * n2)

            # save the x and z components
            @inbounds xz[1] = x
            @inbounds xz[2] = z

            # project velocity onto line of sight
            @inbounds z_rot[i,j] = n2 * angle

            # calculate the limb darkening
            @inbounds ld[i,j] = quad_limb_darkening(μs[i,j], u1, u2)

            # get area element
            @inbounds dA[i,j] = calc_dA_gpu(ρs, ϕc, dϕ, dθ)

            # project onto line of sight
            @inbounds dA[i,j] *= CUDA.abs(a * x + b * y + c * z)
        end
    end
    return nothing
end

function average_subgrid_gpu!(μs_out, μs, wts_out, ld, dA, z_rot_out, z_rot, xz, ax_codes, Nϕ, Nθ)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    # get number of elements along tile side
    k = Int(CUDA.size(μs,1) / Nϕ)
    l = Int(CUDA.size(μs,2) / Nθ)

    # total number of elements output array
    num_tiles = Nϕ *  Nθ

    for t in idx:sdx:num_tiles
        # get index for output array
        row = (t - 1) ÷ Nθ
        col = (t - 1) % Nϕ

        # get indices for input array
        i = row * k + 1
        j = col * l + 1

        # set up sum holders for scalars
        μ_sum = CUDA.zero(CUDA.eltype(μs))
        v_sum = CUDA.zero(CUDA.eltype(z_rot))
        ld_sum = CUDA.zero(CUDA.eltype(ld))
        dA_sum = CUDA.zero(CUDA.eltype(dA))

        # set up sum holders for vector components
        x_sum = CUDA.zero(CUDA.eltype(xz))
        z_sum = CUDA.zero(CUDA.eltype(xz))

        # initiate counter
        count = 0

        for ti in i:i+k-1, tj in j:j+l-1
            if μs[ti, tj] > 0.0
                # sum on scalar quantities
                μ_sum += μs[ti, tj]
                v_sum += z_rot[ti, tj]
                ld_sum += ld[ti, tj]
                dA_sum += dA[ti, tj]

                # sum on vector components
                x_sum += xz[ti, tj, 1]
                z_sum += xz[ti, tj, 2]

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
            @inbounds dat_idx[i,j] = idx
            @inbounds z_cbs[i,j] = cbsall[idx]
        end
    end
    return nothing
end
