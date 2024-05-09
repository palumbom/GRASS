function calc_eclipse_quantities_gpu!(epoch, obs_long, obs_lat, alt, wavelength,
                                        disk::DiskParamsEclipse{T2}, gpu_allocs::GPUAllocsEclipse{T1}) where {T1<:AF, T2<:AF}

    ϕc = gpu_allocs.ϕc
    θc = gpu_allocs.θc
    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot

    # convert scalars from disk params to desired precision
    A = convert(T2, disk.A)
    B = convert(T2, disk.B)
    C = convert(T2, disk.C)
    u1 = convert(T2, disk.u1)
    u2 = convert(T2, disk.u2)

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
    OS_bary_gpu = CuArray(OS_bary)

    # get vector from observatory on earth's surface to moon center
    OM_bary, OM_lt, OM_dlt = spkltc(301, epoch, "J2000", lt_flag, BO_bary)
    OM_bary_gpu = CuArray(OM_bary)

    # get modified epch
    epoch_lt = epoch - OS_lt

    # get rotation matrix for sun
    sun_rot_mat = pxform("IAU_SUN", "J2000", epoch_lt)
    sun_rot_mat_gpu = CuArray(sun_rot_mat)

    #LD - gpu array
    lambda_nm_gpu = CuArray(lambda_nm)
    a0_gpu = CuArray(a0)
    a1_gpu = CuArray(a1)
    a2_gpu = CuArray(a2)
    a3_gpu = CuArray(a3)
    a4_gpu = CuArray(a4)
    a5_gpu = CuArray(a5)

    @cusync begin
        Nθ = CuArray{Float64}(disk.Nθ)
    end

    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(CUDA.length(μs), prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 calc_eclipse_quantities_gpu!(wavelength, ϕc, θc, μs, z_rot,
                                                                                          Nϕ, Nθ, Nsubgrid, Nθ_max, ld, ext, dA,
                                                                                          moon_radius, OS_bary_gpu, OM_bary_gpu, sun_rot_mat_gpu, 
                                                                                          sun_radius, A, B, C, u1, u2,
                                                                                          lambda_nm_gpu, a0_gpu, a1_gpu, a2_gpu, a3_gpu, a4_gpu, a5_gpu)
    return nothing
end                               

function calc_eclipse_quantities_gpu!(wavelength, ϕc, θc, μs, z_rot,
                                        Nϕ, Nθ, Nsubgrid, Nθ_max, ld, ext, dA,
                                        moon_radius, OS_bary, OM_bary, sun_rot_mat, 
                                        sun_radius, A, B, C, u1, u2,
                                        lambda_nm, a0, a1, a2, a3, a4, a5)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    # # get latitude subtile step size
    # N_ϕ_edges = Nϕ * Nsubgrid
    # dϕ = (CUDA.deg2rad(90.0) - CUDA.deg2rad(-90.0)) / (N_ϕ_edges)

    # get number of elements along tile side
    k = Nsubgrid

    # total number of elements output array
    num_tiles = Nϕ * Nθ_max

    # get latitude subtile step size
    N_ϕ_edges = Nϕ * Nsubgrid
    dϕ = π / (N_ϕ_edges)

    # # linear index over course grid tiles
    # for i in idx:sdx:CUDA.length(μs)
    # linear index over course grid tiles
    for t in idx:sdx:num_tiles
        # # get dϕ and dθ for large tiles
        # dϕ_large = (CUDA.deg2rad(90.0) - CUDA.deg2rad(-90.0)) / Nϕ
        # dθ_large = (CUDA.deg2rad(360.0) - CUDA.deg2rad(0.0)) / get_Nθ(ϕc[i], dϕ_large)

        # # get starting edge of large tiles
        # ϕl = ϕc[i] - dϕ_large / 2.0
        # θl = θc[i] - dθ_large / 2.0

        # # get number of longitude tiles in course latitude slice
        # N_θ_edges = get_Nθ(ϕc[i], dϕ_large)
        # dθ = (CUDA.deg2rad(360.0) - CUDA.deg2rad(0.0)) / (N_θ_edges * Nsubgrid)

        # get index for output array - subgridding
        row = (t - 1) ÷ Nθ_max
        col = (t - 1) % Nθ_max
        # @cuprintln("row ", row)
        # @cuprintln("col ", col)
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

        # initiate counter
        μ_count = 0
        count = 0
        
        # # loop over ϕ subgrid
        # for j in 1:Nsubgrid
        #     # get coordinates of latitude subtile center
        #     ϕc_sub = ϕl + (dϕ/2.0) + (j - 1) * dϕ

        #     # loop over θsubgrid
        #     for k in 1:Nsubgrid
        #         # get coordinates of longitude subtile center
        #         θsub = θl + (dθ/2.0) + (k - 1) * dθ

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
                # @cuprintln(typeof(ϕc_sub), " ", typeof(θc_sub))
                # @cuprintln(ϕc_sub, " ", θc_sub)

                # # get cartesian coords of patch
                x, y, z = sphere_to_cart_gpu_eclipse(sun_radius, ϕc_sub, θc_sub)
                #@cuprintln(x," ", y, " ", z)

                # get vector from spherical circle center to surface patch
                a = x
                b = y
                c = CUDA.zero(CUDA.eltype(μs))
                #@cuprintln(a, " ", b, " ", c)

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
                #@cuprintln(rp)

                # set magnitude of vector
                d *= rp
                e *= rp
                f *= rp
                #@cuprintln(d, " ", e, " ", f)

                #xyz rotated for SP bary
                x_new = sun_rot_mat[1] * x + sun_rot_mat[4] * y + sun_rot_mat[7] * z
                y_new = sun_rot_mat[2] * x + sun_rot_mat[5] * y + sun_rot_mat[8] * z
                z_new = sun_rot_mat[3] * x + sun_rot_mat[6] * y + sun_rot_mat[9] * z
                #vel component rotated for SP bary
                vx = sun_rot_mat[1] * d + sun_rot_mat[4] * e + sun_rot_mat[7] * f
                vy = sun_rot_mat[2] * d + sun_rot_mat[5] * e + sun_rot_mat[8] * f
                vz = sun_rot_mat[3] * d + sun_rot_mat[6] * e + sun_rot_mat[9] * f
                #@cuprintln(x_new, " ", y_new, " ", z_new, " ", vx, " ", vy, " ", vz)

                #OP_bary state vector
                OP_bary_x = OS_bary[1] + x_new
                OP_bary_y = OS_bary[2] + y_new
                OP_bary_z = OS_bary[3] + z_new

                OP_bary_vx = OS_bary[4] + vx
                OP_bary_vy = OS_bary[5] + vy
                OP_bary_vz = OS_bary[6] + vz
                #@cuprintln(OP_bary_x, " ", OP_bary_y, " ", OP_bary_z, " ", OP_bary_vx, " ", OP_bary_vy, " ", OP_bary_vz)

                # calculate mu
                μ_sub = calc_mu_gpu(x_new, y_new, z_new, OP_bary_x, OP_bary_y, OP_bary_z) 
                # @cuprintln(μ_sub)
                # @cuprintln(i,j)
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub
                μ_count += 1

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
                d2 = (OM_bary[1] * OP_bary_x + OM_bary[2] * OP_bary_y + OM_bary[3] * OP_bary_z) / (n2 * n1)
                if (d2 < atan(moon_radius/n2))
                    continue
                end

                # get limb darkening
                ld_sub = quad_limb_darkening_gpu_eclipse(μ_sub, wavelength, lambda_nm, a0, a1, a2, a3, a4, a5)
                ld_sum += ld_sub

                z_rot_numerator += z_rot_sub * dA_sub * ld_sub
                z_rot_denominator += dA_sub * ld_sub
            end
        end
        # take averages
        if count > 0
            # set scalar quantity elements as average
            @inbounds μs[i,j] = μ_sum / μ_count
            @cuprintln(i,j)
            @cuprintln(μ_sum / μ_count)
            @inbounds dA[i,j] = dA_sum 
            @inbounds ld[i,j] = ld_sum / μ_count
            @inbounds z_rot[i,j] = z_rot_numerator / z_rot_denominator
        else
            @inbounds μs[i,j] = 0.0
            @inbounds dA[i,j] = 0.0
            @inbounds ld[i,j] = 0.0
            @inbounds z_rot[i,j] = 0.0
        end
    end

    return nothing
end