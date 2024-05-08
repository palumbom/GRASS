function precompute_quantities_gpu!(disk::DiskParams{T1}, μs::CuArray{T2,2}, dA::CuArray{T2,2}, ld::CuArray{T2,3},
                                    wts::CuArray{T2,2}, z_rot::CuArray{T2,2},
                                    ax_codes::CuArray{Int32,2}) where {T1<:AF, T2<:AF}
    # get precision from GPU allocs
    precision = eltype(μs)

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

    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(Nϕ * Nθ_max, prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 precompute_quantities_gpu!(μs, dA, ld, wts, z_rot, ax_codes, Nϕ,
                                                                                       Nθ_max, Nsubgrid, Nθ, R_x, O⃗,
                                                                                       ρs, A, B, C, v0, u1, u2)

    CUDA.synchronize()
    return nothing
end

function precompute_quantities_gpu!(μs, dA, ld, wts, z_rot, ax_codes, Nϕ, Nθ_max, Nsubgrid,
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
    dϕ = π / (N_ϕ_edges)

    # linear index over course grid tiles
    for t in idx:sdx:num_tiles
        # get index for output array
        row = (t - 1) ÷ Nθ_max
        col = (t - 1) % Nθ_max
        m = row + 1
        n = col + 1

        # get indices for input array
        i = row * k + 1
        j = col * k + 1

        # get number of longitude tiles in course latitude slice
        N_θ_edges = Nθ[m] * Nsubgrid

        # set up sum holders for scalars
        μ_sum = CUDA.zero(CUDA.eltype(μs))
        v_sum = CUDA.zero(CUDA.eltype(z_rot))
        ld_sum = CUDA.zero(CUDA.eltype(wts))
        dA_sum = CUDA.zero(CUDA.eltype(wts))

        # set up sum holders for vector components
        x_sum = CUDA.zero(CUDA.eltype(μs))
        y_sum = CUDA.zero(CUDA.eltype(μs))

        # initiate counter
        count = 0

        # loop over latitude sub tiles
        for ti in i:i+k-1
            # get coordinates of latitude subtile center
            ϕc = -π/2 + (dϕ/2.0) + (ti - 1) * dϕ

            # loop over longitude subtiles
            for tj in j:j+k-1
                # move on if we've looped past last longitude
                if tj > N_θ_edges
                    continue
                end

                # get longitude subtile step size
                dθ = 2π / (N_θ_edges)

                # get longitude
                θc = (dθ/2.0) + (tj - 1) * dθ

                # get cartesian coords
                x, y, z = sphere_to_cart_gpu(ρs, ϕc, θc)

                # get vector from spherical circle center to surface patch
                a = x
                b = CUDA.zero(CUDA.eltype(μs))
                c = z

                # take cross product to get vector in direction of rotation
                d = - ρs * c
                e = CUDA.zero(CUDA.eltype(μs))
                f = ρs * a

                # make it a unit vector
                def_norm = CUDA.sqrt(d^2.0 + e^2.0 + f^2.0)
                d /= def_norm
                e /= def_norm
                f /= def_norm

                # set magnitude by differential rotation
                rp = 2π * ρs * CUDA.cos(ϕc) / rotation_period_gpu(ϕc, A, B, C)

                # get in units of c
                rp /= c_Rsun_day

                # set magnitude of vector
                d *= rp
                e *= rp
                f *= rp

                # rotate xyz by inclination
                x, y, z = rotate_vector_gpu(x, y, z, R_x)

                # rotate xyz by inclination and calculate mu
                μ_sub = calc_mu_gpu(x, y, z, O⃗)
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub

                # get limb darkening
                for l in 1:5
                    ld[m,n,l] += quad_limb_darkening_gpu(μ_sub, u1, u2)
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
                v_sum += (n2 * angle)

                # get projected area element
                dA_i = calc_dA_gpu(ρs, ϕc, dϕ, dθ)
                dA_i *= μ_sub
                dA_sum += dA_i

                # sum on vector components
                x_sum += x
                y_sum += y

                # iterate counter
                count += 1
            end
        end

        if count > 0
            # set scalar quantity elements as average
            @inbounds μs[m, n] = μ_sum / count
            @inbounds dA[m, n] = dA_sum
            @inbounds z_rot[m, n] = v_sum / count
            @inbounds wts[m, n] = 0.0

            for l in 1:5
                @inbounds ld[m, n, l] /= count
            end

            # set vector components as average
            @inbounds xx = x_sum / count
            @inbounds yy = y_sum / count

            # get axis code
            @inbounds ax_codes[m, n] = find_nearest_ax_gpu(xx, yy)
        else
            # zero out elements if count is 0 (avoid div by 0)
            @inbounds μs[m, n] = 0.0
            @inbounds dA[m, n] = 0.0
            @inbounds wts[m, n] = 0.0
            @inbounds z_rot[m, n] = 0.0

            for l in 1:5
                @inbounds ld[m, n, l] = 0.0
            end

            # set axis code to 0
            @inbounds ax_codes[m, n] = 0
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

    threads1 = 256
    blocks1 = cld(length(μs), prod(threads1))

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

    for i in idx:sdx:CUDA.length(μs)
        if μs[i] <= 0.0
            continue
        end

        # find the data index for location on disk
        idx = find_data_index_gpu(μs[i], ax_codes[i], disc_mu, disc_ax)
        @inbounds dat_idx[i] = idx
        @inbounds z_cbs[i] = cbsall[idx]
    end
    return nothing
end
