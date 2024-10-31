function calc_eclipse_quantities_gpu!(epoch, obs_long, obs_lat, alt, wavelength,
                                      disk::DiskParamsEclipse{T2}, 
                                      gpu_allocs::GPUAllocsEclipse{T1}) where {T1<:AF, T2<:AF}

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    ax_codes = gpu_allocs.ax_codes

    # re-zero everything
    @cusync begin
        μs .= 0.0
        ld .= 0.0
        ext .= 0.0
        dA .= 0.0
        z_rot .= 0.0
    end

    # convert scalars from disk params to desired precision
    A = convert(T2, disk.A)
    B = convert(T2, disk.B)
    C = convert(T2, disk.C)
    # u1 = convert(T2, disk.u1)
    # u2 = convert(T2, disk.u2)

    # geometry from disk
    #Nϕ = convert(T2, disk.N)
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    Nsubgrid = disk.Nsubgrid

    #query JPL horizons for E, S, M position (km) and velocities (km/s)
    BE_bary = spkssb(399,epoch,"J2000")

    #determine xyz earth coordinates for lat/long of observatory
    flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
    EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)
    #set earth velocity vectors
    EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
    #transform into ICRF frame
    EO_bary = sxform("ITRF93", "J2000", epoch) * EO_earth

    # get vector from barycenter to observatory on Earth's surface
    BO_bary = BE_bary .+ EO_bary

    # set string for ltt and abberation
    lt_flag = "CN+S"

    # get light travel time corrected OS vector
    OS_bary, OS_lt, OS_dlt = spkltc(10, epoch, "J2000", lt_flag, BO_bary)
    @cusync OS_bary_gpu = CuArray(OS_bary)

    # get vector from observatory on earth's surface to moon center
    OM_bary, OM_lt, OM_dlt = spkltc(301, epoch, "J2000", lt_flag, BO_bary)
    @cusync OM_bary_gpu = CuArray(OM_bary)

    # get modified epch
    epoch_lt = epoch - OS_lt

    # get rotation matrix for sun
    sun_rot_mat = pxform("IAU_SUN", "J2000", epoch_lt)
    @cusync sun_rot_mat_gpu = CuArray(sun_rot_mat)

    #LD - gpu array
    lambda_nm_gpu = CuArray(lambda_nm)
    @cusync begin
        a0_gpu = CuArray(a0)
        a1_gpu = CuArray(a1)
        a2_gpu = CuArray(a2)
        a3_gpu = CuArray(a3)
        a4_gpu = CuArray(a4)
        a5_gpu = CuArray(a5)
    end

    @cusync begin
        Nθ = CuArray{Float64}(disk.Nθ)
    end

    @cusync wavelength_gpu = CuArray(wavelength)

    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(CUDA.length(μs), prod(threads1))
    @cusync @cuda threads=threads1 blocks=blocks1 calc_eclipse_quantities_gpu!(wavelength_gpu, μs, z_rot, ax_codes,
                                                                               Nϕ, Nθ, Nsubgrid, Nθ_max, ld, ext, dA,
                                                                               moon_radius, OS_bary_gpu, OM_bary_gpu,
                                                                               sun_rot_mat_gpu, sun_radius, A, B, C,
                                                                               #u1,  u2, 
                                                                               lambda_nm_gpu, a0_gpu, a1_gpu,
                                                                               a2_gpu, a3_gpu, a4_gpu, a5_gpu)
    return nothing
end                               

function calc_eclipse_quantities_gpu!(wavelength, μs, z_rot, ax_codes,
                                        Nϕ, Nθ, Nsubgrid, Nθ_max, ld, ext,
                                        dA, moon_radius, OS_bary, OM_bary,
                                        sun_rot_mat, sun_radius, A, B, C,
                                        #u1, u2, 
                                        lambda_nm, a0, a1, a2, a3,
                                        a4, a5)
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
        z_rot_numerator = CUDA.zero(CUDA.eltype(z_rot))
        z_rot_denominator = CUDA.zero(CUDA.eltype(z_rot))
        ld_sum = CUDA.zero(CUDA.eltype(ld))
        ext_sum = CUDA.zero(CUDA.eltype(ext))
        dA_sum = CUDA.zero(CUDA.eltype(dA))
        x_sum = CUDA.zero(CUDA.eltype(μs))
        y_sum = CUDA.zero(CUDA.eltype(μs))

        # initiate counter
        μ_count = 0
        count = 0

        # loop over latitude sub tiles
        for ti in i:i+k-1
            # get coordinates of latitude subtile center
            ϕc_sub = -π/2 + (dϕ/2.0) + (ti - 1) * dϕ
    
            # loop over longitude subtiles
            for tj in j:j+k-1
                # move on if we've looped past last longitude
                if tj > N_θ_edges
                    continue
                end
    
                # get longitude subtile step size
                dθ = 2π / (N_θ_edges)
    
                # get longitude
                θc_sub = (dθ/2.0) + (tj - 1) * dθ  

                # # get cartesian coords of patch
                x, y, z = sphere_to_cart_gpu_eclipse(sun_radius, ϕc_sub, θc_sub)

                # get vector from spherical circle center to surface patch
                a = x
                b = y
                c = CUDA.zero(CUDA.eltype(μs))

                # take cross product to get vector in direction of rotation
                d = - sun_radius * b
                e = sun_radius * a
                f = CUDA.zero(CUDA.eltype(μs))

                # make it a unit vector
                def_norm = CUDA.sqrt(d^2.0 + e^2.0 + f^2.0)
                d /= def_norm
                e /= def_norm
                f /= def_norm

                # set magnitude by differential rotation
                rp = 2π * sun_radius * CUDA.cos(ϕc_sub) / rotation_period_gpu(ϕc_sub, A, B, C)

                # get in units of c
                rp /= 86400.0

                # set magnitude of vector
                d *= rp
                e *= rp
                f *= rp

                #xyz rotated for SP bary
                x_new = sun_rot_mat[1] * x + sun_rot_mat[4] * y + sun_rot_mat[7] * z
                y_new = sun_rot_mat[2] * x + sun_rot_mat[5] * y + sun_rot_mat[8] * z
                z_new = sun_rot_mat[3] * x + sun_rot_mat[6] * y + sun_rot_mat[9] * z
                #vel component rotated for SP bary
                vx = sun_rot_mat[1] * d + sun_rot_mat[4] * e + sun_rot_mat[7] * f
                vy = sun_rot_mat[2] * d + sun_rot_mat[5] * e + sun_rot_mat[8] * f
                vz = sun_rot_mat[3] * d + sun_rot_mat[6] * e + sun_rot_mat[9] * f

                #OP_bary state vector
                OP_bary_x = OS_bary[1] + x_new
                OP_bary_y = OS_bary[2] + y_new
                OP_bary_z = OS_bary[3] + z_new

                OP_bary_vx = OS_bary[4] + vx
                OP_bary_vy = OS_bary[5] + vy
                OP_bary_vz = OS_bary[6] + vz

                # calculate mu
                μ_sub = calc_mu_gpu(x_new, y_new, z_new, OP_bary_x, OP_bary_y, OP_bary_z) 
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub
                μ_count += 1

                # sum on vector components
                x_sum += x_new
                y_sum += y_new

                # get OP_bary and SP_bary between them and find projected_velocities_no_cb
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n2 = CUDA.sqrt(vx^2.0 + vy^2.0 + vz^2.0)
                angle = (OP_bary_x * vx + OP_bary_y * vy + OP_bary_z * vz) / (n1 * n2)
                v_rot_sub = (n2 * angle)
                v_rot_sub *= 1000.0
                z_rot_sub = v_rot_sub / c_ms

                n2 = CUDA.sqrt(OS_bary[4]^2.0 + OS_bary[5]^2.0 + OS_bary[6]^2.0)
                angle = (OP_bary_x * OS_bary[4] + OP_bary_y * OS_bary[5] + OP_bary_z * OS_bary[6]) / (n1 * n2)
                v_orbit_sub = (n2 * angle)
                v_orbit_sub *= 1000.0
                z_rot_sub += (v_orbit_sub / c_ms)

                # get projected area element
                dA_sub = calc_dA_gpu(sun_radius, ϕc_sub, dϕ, dθ)
                dA_sub *= μ_sub
                dA_sum += dA_sub

                # iterate counter 
                count += 1

                #calculate distance
                n2 = CUDA.sqrt(OM_bary[1]^2.0 + OM_bary[2]^2.0 + OM_bary[3]^2.0)  
                d2 = acos((OM_bary[1] * OP_bary_x + OM_bary[2] * OP_bary_y + OM_bary[3] * OP_bary_z) / (n2 * n1))
                if (d2 < atan(moon_radius/n2))
                    continue
                end

                for wl in eachindex(wavelength)
                    # get limb darkening
                    ld[m,n,wl] += quad_limb_darkening_gpu_eclipse(μ_sub, wavelength[wl], lambda_nm, a0, a1, a2, a3, a4, a5)
                end
                z_rot_numerator += z_rot_sub * dA_sub * ld[m,n,1]
                z_rot_denominator += dA_sub * ld[m,n,1]
            end
        end
        # take averages
        if count > 0
            # set scalar quantity elements as average
            @inbounds μs[m,n] = μ_sum / μ_count
            @inbounds dA[m,n] = dA_sum 
            for wl in eachindex(wavelength)
                @inbounds ld[m,n,wl] /= count
            end
            @inbounds z_rot[m,n] = z_rot_numerator / z_rot_denominator

            # set vector components as average
            @inbounds xx = x_sum / μ_count
            @inbounds yy = y_sum / μ_count

            @inbounds ax_codes[m, n] = find_nearest_ax_gpu(xx / sun_radius, yy / sun_radius)
        else
            @inbounds μs[m,n] = 0.0
            @inbounds dA[m,n] = 0.0
            for wl in eachindex(wavelength)
                @inbounds ld[m,n,wl] = 0.0
            end
            @inbounds z_rot[m,n] = 0.0
        end
    end

    return nothing
end