function precompute_quantities_gpu!(disk::DiskParams{T1}, gpu_allocs::GPUAllocs{T2}) where {T1<:AF, T2<:AF}
    # get precision from GPU allocs
    precision = eltype(gpu_allocs.λs)

    # convert scalars from disk params to desired precision
    ρs = convert(precision, disk.ρs)
    A = convert(precision, disk.A)
    B = convert(precision, disk.B)
    C = convert(precision, disk.C)
    v0 = convert(precision, disk.v0)
    u1 = convert(precision, disk.u1)
    u2 = convert(precision, disk.u2)

    # get size of sub-tiled grid
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    Nsubgrid = disk.Nsubgrid
    Nϕ_sub = Nϕ * Nsubgrid
    Nθ_sub = maximum(disk.Nθ) * Nsubgrid

    # copy data to GPU
    @cusync begin
        # get observer vectoir and rotation matrix
        O⃗ = CuArray{precision}(disk.O⃗)
        Nθ = CuArray{Int32}(disk.Nθ)
        R_x = CuArray{precision}(disk.R_x)
    end

    # alias from GPU allocs
    μs_out = gpu_allocs.μs
    wts_out = gpu_allocs.wts
    z_rot_out = gpu_allocs.z_rot
    ax_codes = gpu_allocs.ax_codes

    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(Nϕ * Nθ_max, prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 precompute_quantities_gpu!(μs_out, wts_out, z_rot_out,
                                                                                       ax_codes, Nϕ, Nθ_max, Nsubgrid,
                                                                                       Nθ, R_x, O⃗, ρs, A, B, C, v0, u1, u2)

    CUDA.synchronize()
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

function precompute_quantities_gpu!(μs, wts, z_rot, ax_codes, Nϕ, Nθ_max, Nsubgrid,
                                    Nθ, R_x, O⃗, ρs, A, B, C, v0, u1, u2)
# get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    # get number of elements along tile side
    k = Nsubgrid

    # total number of elements output array
    num_tiles = Nϕ * Nθ_max

    # get latitude subtile step size
    N_ϕ_edges = Nϕ * Nsubgrid
    dϕ = (CUDA.deg2rad(90.0) - CUDA.deg2rad(-90.0)) / (N_ϕ_edges)

    # linear index over course grid tiles
    for t in idx:sdx:num_tiles
        # get index for output array
        row = (t - 1) ÷ Nθ_max
        col = (t - 1) % Nθ_max

        # get indices for input array
        i = row * k + 1
        j = col * k + 1

        # set up sum holders for scalars
        μ_sum = CUDA.zero(CUDA.eltype(μs))
        v_sum = CUDA.zero(CUDA.eltype(z_rot))
        ld_sum = CUDA.zero(CUDA.eltype(wts))
        dA_sum = CUDA.zero(CUDA.eltype(wts))

        # set up sum holders for vector components
        x_sum = CUDA.zero(CUDA.eltype(μs))
        z_sum = CUDA.zero(CUDA.eltype(μs))

        # initiate counter
        count = 0

        # get number of longitude tiles in course latitude slice
        m = row + 1
        N_θ_edges = Nθ[m] * Nsubgrid

        # loop over latitude sub tiles
        for ti in i:i+k-1
            # get coordinates of latitude subtile center
            ϕc = CUDA.deg2rad(-90.0) + (dϕ/2.0) + (ti - 1) * dϕ

            # loop over longitude subtiles
            for tj in j:j+k-1
                # move on if we've looped past last longitude
                if tj > N_θ_edges
                    continue
                end

                # get longitude subtile step size
                dθ = (CUDA.deg2rad(360.0) - CUDA.deg2rad(0.0)) / (N_θ_edges)

                # get longitude
                θc = CUDA.deg2rad(0.0) + (dθ/2.0) + (tj - 1) * dθ

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
                rp = (v0 / rotation_period_gpu(ϕc, A, B, C))
                # rp = 33950.0/3e8 * (1.0 - A * sin(ϕc)^2.0)
                d *= rp
                e *= rp

                # rotate xyz by inclination
                x, y, z = rotate_vector_gpu(x, y, z, R_x)

                # rotate xyz by inclination and calculate mu
                μ_sub = calc_mu_gpu(x, y, z, O⃗)
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub

                # get limb darkening
                ld_sum += quad_limb_darkening_gpu(μ_sub, u1, u2)

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
                v_sum += (n2 * angle)

                # get projected area element
                dA = calc_dA_gpu(ρs, ϕc, dϕ, dθ)
                dA *= CUDA.abs(a * x + b * y + c * z)
                dA *= CUDA.sqrt(O⃗[1]^2.0 + O⃗[2]^2.0 + O⃗[3]^2.0)
                dA_sum += dA

                # sum on vector components
                x_sum += x
                z_sum += z

                # iterate counter
                count += 1
            end
        end

        if count > 0
            # set scalar quantity elements as average
            @inbounds μs[row + 1, col + 1] = μ_sum / count
            @inbounds z_rot[row + 1, col + 1] = v_sum / count
            @inbounds wts[row + 1, col + 1] = (ld_sum / count) * dA_sum

            # set scalar quantity elements as average
            @inbounds xx = x_sum / count
            @inbounds zz = z_sum / count

            # get axis code
            @inbounds ax_codes[row + 1, col + 1] = find_nearest_ax_gpu(xx, zz)
        else
            # zero out elements if count is 0 (avoid div by 0)
            @inbounds μs[row + 1, col + 1] = 0.0
            @inbounds wts[row + 1, col + 1] = 0.0
            @inbounds z_rot[row + 1, col + 1] = 0.0

            # set axis code to 0
            @inbounds ax_codes[row + 1, col + 1] = 0
        end
    end

end
