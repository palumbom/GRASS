function calc_eclipse_quantities_gpu!(epoch::T1, obs_long::T1, obs_lat::T1, alt::T1, wavelength, 
                                      LD_law::String, ext_toggle::T1, ext_coeff::T1,
                                      disk::DiskParamsEclipse{T2}, gpu_allocs::GPUAllocsEclipse{T1}) where {T1<:AF, T2<:AF}

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    projected_v = gpu_allocs.projected_v
    earth_v = gpu_allocs.earth_v
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    ax_codes = gpu_allocs.ax_codes

    # re-zero everything
    @cusync begin
        μs .= 0.0
        ld .= 0.0
        projected_v .= 0.0
        earth_v .= 0.0
        ext .= 0.0
        dA .= 0.0
        z_rot .= 0.0
    end

    # convert scalars from disk params to desired precision
    A = convert(T2, disk.A)
    B = convert(T2, disk.B)
    C = convert(T2, disk.C)
    ext_coeff_gpu = convert(T2, ext_coeff)

    if LD_law == "300"
        filtered_df = quad_ld_coeff_300[quad_ld_coeff_300.wavelength .== wavelength, :]
        u1 = convert(T2, filtered_df.law_4u1[1])
        u2 = convert(T2, filtered_df.law_4u2[1])
    end

    if LD_law == "HD"
        filtered_df = quad_ld_coeff_HD[quad_ld_coeff_HD.wavelength .== wavelength, :]
        u1 = convert(T2, filtered_df.law_4u1[1])
        u2 = convert(T2, filtered_df.law_4u2[1])
    end

    if LD_law == "SSD_quadratic"
        filtered_df = quad_ld_coeff_SSD[quad_ld_coeff_SSD.wavelength .== wavelength, :]
        u1 = convert(T2, filtered_df.u1[1])
        u2 = convert(T2, filtered_df.u2[1])
        u3 = NaN
        u4 = NaN
    end

    if LD_law == "SSD_4parameter"
        filtered_df = quad_ld_coeff_SSD[quad_ld_coeff_SSD.wavelength .== wavelength, :]
        u1 = convert(T2, filtered_df.law_4u1[1])
        u2 = convert(T2, filtered_df.law_4u2[1])
        u3 = convert(T2, filtered_df.law_4u3[1])
        u4 = convert(T2, filtered_df.law_4u4[1])
    end

    # geometry from disk
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    Nsubgrid = disk.Nsubgrid

    # query SPICE for E, S, M position (km) and velocities (km/s)
    BE_bary = spkssb(399,epoch,"J2000")

    # determine xyz earth coordinates for lat/long of observatory
    flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
    EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)
    # set earth velocity vectors
    EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
    # transform into ICRF frame
    EO_bary = sxform("IAU_EARTH", "J2000", epoch) * EO_earth
    @cusync EO_bary_gpu = CuArray(EO_bary)

    # get vector from barycenter to observatory on Earth's surface
    BO_bary = BE_bary .+ EO_bary

    # set string for ltt and abberation
    lt_flag = "CN+S"

    @cusync moon_radius_gpu = CuArray([moon_radius])
    @cusync sun_radius_gpu = CuArray([sun_radius])

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

    @cusync begin
        Nθ = CuArray{Float64}(disk.Nθ)
    end

    @cusync wavelength_gpu = CuArray(wavelength)
    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(prod(CUDA.size(μs)), prod(threads1))
    @cusync @cuda threads=threads1 blocks=blocks1 calc_eclipse_quantities_gpu!(wavelength_gpu, μs, z_rot, ax_codes,
                                                                               Nϕ, Nθ, Nsubgrid, Nθ_max, ld, projected_v, earth_v, ext, 
                                                                               dA, moon_radius_gpu, OS_bary_gpu, OM_bary_gpu, EO_bary_gpu,
                                                                               sun_rot_mat_gpu, sun_radius_gpu, A, B, C, u1, u2, u3, u4,
                                                                               ext_toggle, ext_coeff_gpu)
    return
end                               

function calc_eclipse_quantities_gpu!(wavelength, μs, z_rot, ax_codes,
                                      Nϕ, Nθ, Nsubgrid, Nθ_max, ld, projected_v, earth_v, ext,
                                      dA, moon_radius, OS_bary, OM_bary, EO_bary,
                                      sun_rot_mat, sun_radius, A, B, C, u1, u2, u3, u4,
                                      ext_toggle, ext_coeff_gpu)
    sun_radius = sun_radius[1] 
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
        z_rot_numerator = CUDA.zero(CUDA.eltype(μs))
        z_rot_denominator = CUDA.zero(CUDA.eltype(μs))
        dA_sum = CUDA.zero(CUDA.eltype(μs))
        ld_sum = CUDA.zero(CUDA.eltype(μs))
        projected_v_sum = CUDA.zero(CUDA.eltype(μs))
        earth_v_sum = CUDA.zero(CUDA.eltype(μs))
        ext_sum = CUDA.zero(CUDA.eltype(μs))
        x_sum = CUDA.zero(CUDA.eltype(μs))
        y_sum = CUDA.zero(CUDA.eltype(μs))
        z_sum = CUDA.zero(CUDA.eltype(μs))

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

                # get cartesian coords of patch
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

                # xyz rotated for SP bary
                x_new = sun_rot_mat[1] * x + sun_rot_mat[4] * y + sun_rot_mat[7] * z
                y_new = sun_rot_mat[2] * x + sun_rot_mat[5] * y + sun_rot_mat[8] * z
                z_new = sun_rot_mat[3] * x + sun_rot_mat[6] * y + sun_rot_mat[9] * z

                # vel component rotated for SP bary
                vx = sun_rot_mat[1] * d + sun_rot_mat[4] * e + sun_rot_mat[7] * f
                vy = sun_rot_mat[2] * d + sun_rot_mat[5] * e + sun_rot_mat[8] * f
                vz = sun_rot_mat[3] * d + sun_rot_mat[6] * e + sun_rot_mat[9] * f

                # OP_bary state vector
                OP_bary_x = OS_bary[1] + x_new
                OP_bary_y = OS_bary[2] + y_new
                OP_bary_z = OS_bary[3] + z_new

                # calculate mu
                μ_sub = calc_mu_gpu(x_new, y_new, z_new, OP_bary_x, OP_bary_y, OP_bary_z) 
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub
                μ_count += 1

                # sum on vector components
                x_sum += x
                y_sum += y
                z_sum += z

                # get OP_bary and SP_bary between them and find projected_velocities_no_cb
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n2 = CUDA.sqrt(vx^2.0 + vy^2.0 + vz^2.0)
                angle = (OP_bary_x * vx + OP_bary_y * vy + OP_bary_z * vz) / (n1 * n2)
                v_rot_sub = (n2 * angle)
                v_rot_sub *= 1000.0

                n2 = CUDA.sqrt(OS_bary[4]^2.0 + OS_bary[5]^2.0 + OS_bary[6]^2.0)
                angle = (OP_bary_x * OS_bary[4] + OP_bary_y * OS_bary[5] + OP_bary_z * OS_bary[6]) / (n1 * n2)
                v_orbit_sub = (n2 * angle)
                v_orbit_sub *= 1000.0

                # get projected area element
                dA_sub = calc_dA_gpu(sun_radius, ϕc_sub, dϕ, dθ)
                dA_sub *= μ_sub
                dA_sum += dA_sub

                # calculate distance (Solar Eclipse)
                n2 = CUDA.sqrt(OM_bary[1]^2.0 + OM_bary[2]^2.0 + OM_bary[3]^2.0)  
                d2 = acos((OM_bary[1] * OP_bary_x + OM_bary[2] * OP_bary_y + OM_bary[3] * OP_bary_z) / (n2 * n1))
                if (d2 < atan(moon_radius[1]/n2))
                    continue
                end

                z_rot_sub = v_rot_sub / c_ms
                z_rot_sub += (v_orbit_sub / c_ms)

                projected_v_sum += v_rot_sub
                earth_v_sum += v_orbit_sub

                # iterate counter 
                count += 1

                # zenith
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n3 = CUDA.sqrt(EO_bary[1]^2.0 + EO_bary[2]^2.0 + EO_bary[3]^2.0) 
                zenith = (CUDA.acos((EO_bary[1] * OP_bary_x + EO_bary[2] * OP_bary_y + EO_bary[3] * OP_bary_z) / (n3 * n1)))

                for wl in eachindex(wavelength)
                    # get limb darkening
                    if isnan(u3)
                        ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2)
                    else
                        ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2, u3, u4)
                    end
                    ld_sum += ld_sub

                    if CUDA.isone(ext_toggle)
                        ext_sub = CUDA.exp(-((1/(CUDA.cos(zenith)))*ext_coeff_gpu))
                        ext_sum += ext_sub
                    end
                end

                # recalculate value since local vars from wl for loop are lost
                if isnan(u3)
                    ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2)
                else
                    ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2, u3, u4)
                end
                ext_sub = CUDA.exp(-((1/(CUDA.cos(zenith)))*ext_coeff_gpu))

                if ext_toggle == 1.0
                    z_rot_numerator += z_rot_sub * dA_sub * ld_sub * ext_sub
                    z_rot_denominator += dA_sub * ld_sub * ext_sub
                else
                    z_rot_numerator += z_rot_sub * dA_sub * ld_sub
                    z_rot_denominator += dA_sub * ld_sub 
                end
            end
        end
        # take averages
        @inbounds μs[m,n] = μ_sum / μ_count
        @inbounds dA[m,n] = dA_sum 

        if count > 0
            # set scalar quantity elements as average
            @inbounds projected_v[m,n] = projected_v_sum / count
            @inbounds earth_v[m,n] = earth_v_sum / count

            @inbounds z_rot[m,n] = z_rot_numerator / z_rot_denominator

            for wl in eachindex(wavelength)
                @inbounds ld[m,n,wl] = ld_sum / count
                @inbounds ext[m,n,wl] = ext_sum / count
            end

            # set vector components as average
            @inbounds xx = x_sum / μ_count
            @inbounds yy = y_sum / μ_count
            @inbounds zz = z_sum / μ_count

            @inbounds ax_codes[m, n] = find_nearest_ax_gpu(yy / sun_radius, zz / sun_radius)
        end
    end

    return nothing
end

function calc_eclipse_quantities_gpu!(epoch::T1, t, obs_long::T1, obs_lat::T1, alt::T1, wavelength, 
                                      ext_toggle::T1, ext_coeff::T1,
                                      disk::DiskParamsEclipse{T2}, gpu_allocs::GPUAllocsEclipse{T1}, LD1, LD2, LD3, LD4, CB1, CB2, CB3, MF1, MF2) where {T1<:AF, T2<:AF}
    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    projected_v = gpu_allocs.projected_v
    earth_v = gpu_allocs.earth_v
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    ax_codes = gpu_allocs.ax_codes

    # re-zero everything
    @cusync begin
        μs .= 0.0
        ld .= 0.0
        projected_v .= 0.0
        earth_v .= 0.0
        ext .= 0.0
        dA .= 0.0
        z_rot .= 0.0
    end

    # convert scalars from disk params to desired precision
    A = convert(T2, disk.A)
    B = convert(T2, disk.B)
    C = convert(T2, disk.C)
    ext_coeff_gpu = convert(T2, ext_coeff)

    # u1 = convert(T2, LD1)
    # u2 = convert(T2, LD2)
    # u3 = convert(T2, LD3)
    # u4 = convert(T2, LD4)
    filtered_df = quad_ld_coeff_SSD[quad_ld_coeff_SSD.wavelength .== wavelength, :]
    u1 = convert(T2, filtered_df.law_4u1[1])
    u2 = convert(T2, filtered_df.law_4u2[1])
    u3 = convert(T2, filtered_df.law_4u3[1])
    u4 = convert(T2, filtered_df.law_4u4[1])

    CB1 = convert(T2, CB1)
    CB2 = convert(T2, CB2)
    CB3 = convert(T2, CB3)

    MF1 = convert(T2, MF1)
    MF2 = convert(T2, MF2)

    # geometry from disk
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    Nsubgrid = disk.Nsubgrid

    # query SPICE for E, S, M position (km) and velocities (km/s)
    BE_bary = spkssb(399,epoch,"J2000")

    # determine xyz earth coordinates for lat/long of observatory
    flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
    EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)
    # set earth velocity vectors
    EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
    # transform into ICRF frame
    EO_bary = sxform("IAU_EARTH", "J2000", epoch) * EO_earth
    @cusync EO_bary_gpu = CuArray(EO_bary)

    # get vector from barycenter to observatory on Earth's surface
    BO_bary = BE_bary .+ EO_bary

    # set string for ltt and abberation
    lt_flag = "CN+S"

    @cusync moon_radius_gpu = CuArray([moon_radius])
    @cusync earth_radius_gpu = CuArray([earth_radius])
    @cusync sun_radius_gpu = CuArray([sun_radius])

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

    @cusync begin
        Nθ = CuArray{Float64}(disk.Nθ)
    end

    # spots_ra_gpu = CuArray(spots_ra_time_arr[t])
    # spots_dec_gpu = CuArray(spots_dec_time_arr[t])
    spot_lon = CuArray(deg2rad.(spots_info.lon))
    spot_lat = CuArray(deg2rad.(spots_info.lat))

    spot_xyz = map((x,y) -> sphere_to_cart_eclipse.(sun_radius, x, y), deg2rad.(spots_info.lat), deg2rad.(spots_info.lon))
    spot_xyz = map(x -> (sun_rot_mat * x) .+ OS_bary[1:3], spot_xyz)
    spot_xyz_mat = reduce(hcat, spot_xyz)
    spot_xyz_gpu = CuArray(spot_xyz_mat)
    
    contrast = CuArray(spots_info.contrast)
    diameter_arcsec = CuArray(spots_info.diameter_arcsec)
    diameter_km = CuArray(spots_info.diameter_km)

    @cusync wavelength_gpu = CuArray(wavelength)
    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(prod(CUDA.size(μs)), prod(threads1))
    @cusync @cuda threads=threads1 blocks=blocks1 calc_eclipse_quantities_gpu!(wavelength_gpu, μs, z_rot, ax_codes,
                                                                               Nϕ, Nθ, Nsubgrid, Nθ_max, ld, projected_v, earth_v, ext, 
                                                                               dA, moon_radius_gpu, OS_bary_gpu, OM_bary_gpu, EO_bary_gpu,
                                                                               sun_rot_mat_gpu, sun_radius_gpu, A, B, C, u1, u2, u3, u4, CB1, CB2, CB3, MF1, MF2,
                                                                               ext_toggle, ext_coeff_gpu, earth_radius_gpu, spot_xyz_gpu, contrast, diameter_arcsec, diameter_km)
    return
end                               

function calc_eclipse_quantities_gpu!(wavelength, μs, z_rot, ax_codes,
                                      Nϕ, Nθ, Nsubgrid, Nθ_max, ld, projected_v, earth_v, ext,
                                      dA, moon_radius, OS_bary, OM_bary, EO_bary,
                                      sun_rot_mat, sun_radius, A, B, C, u1, u2, u3, u4, CB1, CB2, CB3, MF1, MF2,
                                      ext_toggle, ext_coeff_gpu, earth_radius, spot_xyz_gpu, contrast, diameter_arcsec, diameter_km)
    sun_radius = sun_radius[1] 
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
        z_rot_numerator = CUDA.zero(CUDA.eltype(μs))
        z_rot_denominator = CUDA.zero(CUDA.eltype(μs))
        dA_sum = CUDA.zero(CUDA.eltype(μs))
        ld_sum = CUDA.zero(CUDA.eltype(μs))
        projected_v_sum = CUDA.zero(CUDA.eltype(μs))
        earth_v_sum = CUDA.zero(CUDA.eltype(μs))
        ext_sum = CUDA.zero(CUDA.eltype(μs))
        x_sum = CUDA.zero(CUDA.eltype(μs))
        y_sum = CUDA.zero(CUDA.eltype(μs))
        z_sum = CUDA.zero(CUDA.eltype(μs))

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

                # get cartesian coords of patch
                x, y, z = sphere_to_cart_gpu_eclipse(sun_radius, ϕc_sub, θc_sub)

                # spot1_x, spot1_y, spot1_z = sphere_to_cart_gpu_eclipse(sun_radius, spot_lat[1], spot_lon[1])
                # spot2_x, spot2_y, spot2_z = sphere_to_cart_gpu_eclipse(sun_radius, spot_lat[2], spot_lon[2])
                # spot3_x, spot3_y, spot3_z = sphere_to_cart_gpu_eclipse(sun_radius, spot_lat[3], spot_lon[3])

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

                # xyz rotated for SP bary
                x_new = sun_rot_mat[1] * x + sun_rot_mat[4] * y + sun_rot_mat[7] * z
                y_new = sun_rot_mat[2] * x + sun_rot_mat[5] * y + sun_rot_mat[8] * z
                z_new = sun_rot_mat[3] * x + sun_rot_mat[6] * y + sun_rot_mat[9] * z

                # spot1_x_new = sun_rot_mat[1] * spot1_x + sun_rot_mat[4] * spot1_y + sun_rot_mat[7] * spot1_z
                # spot1_y_new = sun_rot_mat[2] * spot1_x + sun_rot_mat[5] * spot1_y + sun_rot_mat[8] * spot1_z
                # spot1_z_new = sun_rot_mat[3] * spot1_x + sun_rot_mat[6] * spot1_y + sun_rot_mat[9] * spot1_z

                # spot2_x_new = sun_rot_mat[1] * spot2_x + sun_rot_mat[4] * spot2_y + sun_rot_mat[7] * spot2_z
                # spot2_y_new = sun_rot_mat[2] * spot2_x + sun_rot_mat[5] * spot2_y + sun_rot_mat[8] * spot2_z
                # spot2_z_new = sun_rot_mat[3] * spot2_x + sun_rot_mat[6] * spot2_y + sun_rot_mat[9] * spot2_z

                # spot3_x_new = sun_rot_mat[1] * spot3_x + sun_rot_mat[4] * spot3_y + sun_rot_mat[7] * spot3_z
                # spot3_y_new = sun_rot_mat[2] * spot3_x + sun_rot_mat[5] * spot3_y + sun_rot_mat[8] * spot3_z
                # spot3_z_new = sun_rot_mat[3] * spot3_x + sun_rot_mat[6] * spot3_y + sun_rot_mat[9] * spot3_z

                # vel component rotated for SP bary
                vx = sun_rot_mat[1] * d + sun_rot_mat[4] * e + sun_rot_mat[7] * f
                vy = sun_rot_mat[2] * d + sun_rot_mat[5] * e + sun_rot_mat[8] * f
                vz = sun_rot_mat[3] * d + sun_rot_mat[6] * e + sun_rot_mat[9] * f

                # OP_bary state vector
                OP_bary_x = OS_bary[1] + x_new
                OP_bary_y = OS_bary[2] + y_new
                OP_bary_z = OS_bary[3] + z_new

                # Ospot1_bary_x = OS_bary[1] + spot1_x_new
                # Ospot1_bary_y = OS_bary[2] + spot1_y_new
                # Ospot1_bary_z = OS_bary[3] + spot1_z_new

                # Ospot2_bary_x = OS_bary[1] + spot2_x_new
                # Ospot2_bary_y = OS_bary[2] + spot2_y_new
                # Ospot2_bary_z = OS_bary[3] + spot2_z_new

                # Ospot3_bary_x = OS_bary[1] + spot3_x_new
                # Ospot3_bary_y = OS_bary[2] + spot3_y_new
                # Ospot3_bary_z = OS_bary[3] + spot3_z_new

                # calculate mu
                μ_sub = calc_mu_gpu(x_new, y_new, z_new, OP_bary_x, OP_bary_y, OP_bary_z) 
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub
                μ_count += 1

                # sum on vector components
                x_sum += x
                y_sum += y
                z_sum += z

                # get OP_bary and SP_bary between them and find projected_velocities_no_cb
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n2 = CUDA.sqrt(vx^2.0 + vy^2.0 + vz^2.0)
                angle = (OP_bary_x * vx + OP_bary_y * vy + OP_bary_z * vz) / (n1 * n2)
                v_rot_sub = (n2 * angle) + (CB1*CUDA.exp(μ_sub)^2.0 + CB2*CUDA.exp(μ_sub) + CB3) + MF1*meridional_terms(rad2deg(ϕc_sub), rad2deg(θc_sub), 5.94)[1] + MF2*meridional_terms(rad2deg(ϕc_sub), rad2deg(θc_sub), 5.94)[2]
                v_rot_sub_no_cb = (n2 * angle)
                v_rot_sub *= 1000.0
                v_rot_sub_no_cb *= 1000.0

                n2 = CUDA.sqrt(OS_bary[4]^2.0 + OS_bary[5]^2.0 + OS_bary[6]^2.0)
                angle = (OP_bary_x * OS_bary[4] + OP_bary_y * OS_bary[5] + OP_bary_z * OS_bary[6]) / (n1 * n2)
                v_orbit_sub = (n2 * angle)
                v_orbit_sub *= 1000.0

                # get projected area element
                dA_sub = calc_dA_gpu(sun_radius, ϕc_sub, dϕ, dθ)
                dA_sub *= μ_sub
                dA_sum += dA_sub

                # calculate distance (Solar Eclipse)
                n2 = CUDA.sqrt(OM_bary[1]^2.0 + OM_bary[2]^2.0 + OM_bary[3]^2.0)  
                d2 = acos((OM_bary[1] * OP_bary_x + OM_bary[2] * OP_bary_y + OM_bary[3] * OP_bary_z) / (n2 * n1))
                if (d2 < atan(moon_radius[1]/n2))
                    continue
                end

                # iterate counter 
                count += 1

                z_rot_sub = v_rot_sub / c_ms
                z_rot_sub += (v_orbit_sub / c_ms)

                projected_v_sum += v_rot_sub
                earth_v_sum += v_orbit_sub

                # zenith
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n3 = CUDA.sqrt(EO_bary[1]^2.0 + EO_bary[2]^2.0 + EO_bary[3]^2.0) 
                zenith = (CUDA.acos((EO_bary[1] * OP_bary_x + EO_bary[2] * OP_bary_y + EO_bary[3] * OP_bary_z) / (n3 * n1)))

                # skip_patch = false
                # contrast_value = 0.0
                # for s in 1:size(spot_xyz_gpu, 2)
                #     n2 = CUDA.sqrt(spot_xyz_gpu[1, s]^2.0 + spot_xyz_gpu[2, s]^2.0 + spot_xyz_gpu[3, s]^2.0)  
                #     d2 = acos((spot_xyz_gpu[1, s] * OP_bary_x + spot_xyz_gpu[2, s] * OP_bary_y + spot_xyz_gpu[3, s] * OP_bary_z) / (n2 * n1))
                #     if (d2 < atan((diameter_km[s]/2)/n2))
                #         skip_patch = true
                #         contrast_value = contrast[s]
                #         break 
                #     end
                # end

                # # skip_patch = false
                # # contrast_value = 0.0
                # # for s in 1:length(contrast)
                # #     n2 = CUDA.sqrt(Ospot1_bary_x^2.0 + Ospot1_bary_y^2.0 + Ospot1_bary_z^2.0)  
                # #     d2 = acos((Ospot1_bary_x * OP_bary_x + Ospot1_bary_y * OP_bary_y + Ospot1_bary_z * OP_bary_z) / (n2 * n1))
                # #     if (d2 < atan((diameter_km[s]/2)/n2))
                # #         skip_patch = true
                # #         contrast_value = contrast[s]
                # #         break 
                # #     end

                # #     n2 = CUDA.sqrt(Ospot2_bary_x^2.0 + Ospot2_bary_y^2.0 + Ospot2_bary_z^2.0)  
                # #     d2 = acos((Ospot2_bary_x * OP_bary_x + Ospot2_bary_y * OP_bary_y + Ospot2_bary_z * OP_bary_z) / (n2 * n1))
                # #     if (d2 < atan((diameter_km[s]/2)/n2))
                # #         skip_patch = true
                # #         contrast_value = contrast[s]
                # #         break 
                # #     end

                # #     n2 = CUDA.sqrt(Ospot3_bary_x^2.0 + Ospot3_bary_y^2.0 + Ospot3_bary_z^2.0)  
                # #     d2 = acos((Ospot3_bary_x * OP_bary_x + Ospot3_bary_y * OP_bary_y + Ospot3_bary_z * OP_bary_z) / (n2 * n1))
                # #     if (d2 < atan((diameter_km[s]/2)/n2))
                # #         skip_patch = true
                # #         contrast_value = contrast[s]
                # #         break 
                # #     end
                # # end

                # if !skip_patch
                #     projected_v_sum += v_rot_sub
                # end
                # if skip_patch
                #     projected_v_sum += v_rot_sub_no_cb
                # end

                for wl in eachindex(wavelength)
                    # get limb darkening
                    ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2, u3, u4)
                    ld_sum += ld_sub
                    # if !skip_patch
                    #     ld_sum += ld_sub
                    # end
                    # if skip_patch
                    #     ld_sum += contrast_value*ld_sub
                    # end

                    if CUDA.isone(ext_toggle)
                        ext_sub = CUDA.exp(-((1/(CUDA.cos(zenith)))*ext_coeff_gpu))
                        ext_sum += ext_sub
                    end
                end

                # recalculate value since local vars from wl for loop are lost
                # if !skip_patch
                #     ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2, u3, u4)
                # end
                # if skip_patch
                #     ld_sub = contrast_value*quad_limb_darkening_gpu(μ_sub, u1, u2, u3, u4)
                # end
                ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2, u3, u4)
                ext_sub = CUDA.exp(-((1/(CUDA.cos(zenith)))*ext_coeff_gpu))

                if ext_toggle == 1.0
                    z_rot_numerator += z_rot_sub * dA_sub * ld_sub * ext_sub
                    z_rot_denominator += dA_sub * ld_sub * ext_sub
                else
                    z_rot_numerator += z_rot_sub * dA_sub * ld_sub
                    z_rot_denominator += dA_sub * ld_sub 
                end
            end
        end
        # take averages
        @inbounds μs[m,n] = μ_sum / μ_count
        @inbounds dA[m,n] = dA_sum 

        if count > 0
            # set scalar quantity elements as average
            @inbounds projected_v[m,n] = projected_v_sum / count
            @inbounds earth_v[m,n] = earth_v_sum / count

            @inbounds z_rot[m,n] = z_rot_numerator / z_rot_denominator

            for wl in eachindex(wavelength)
                @inbounds ld[m,n,wl] = ld_sum / count
                @inbounds ext[m,n,wl] = ext_sum / count
            end

            # set vector components as average
            @inbounds xx = x_sum / μ_count
            @inbounds yy = y_sum / μ_count
            @inbounds zz = z_sum / μ_count

            @inbounds ax_codes[m, n] = find_nearest_ax_gpu(yy / sun_radius, zz / sun_radius)
        end
    end

    return nothing
end

function calc_grid_edge_xyz(epoch, obs_long, obs_lat, alt,
                            disk::DiskParamsEclipse{T2}, 
                            gpu_allocs::GPUAllocsEclipse{T1}) where {T1<:AF, T2<:AF}
    # get phi edges
    ϕe = repeat(collect(disk.ϕe), inner=2)[2:end-1]
    ϕe = CuArray(ϕe)

    # get number of theta edges
    Nθ = repeat(disk.Nθ, inner=2) .+ 1
    θe = zeros(length(ϕe), maximum(Nθ))
    for i in eachindex(ϕe)
        θe[i, 1:Nθ[i]] .= range(0.0, 2π, length=Nθ[i])
    end
    θe = CuArray(θe)

    # allocate memory for cartesian coords
    xs = CUDA.zeros(Float64, CUDA.length(ϕe), CUDA.size(θe,2))
    ys = CUDA.zeros(Float64, CUDA.length(ϕe), CUDA.size(θe,2))
    zs = CUDA.zeros(Float64, CUDA.length(ϕe), CUDA.size(θe,2))

    # query JPL horizons for E, S, M position (km) and velocities (km/s)
    BE_bary = spkssb(399,epoch,"J2000")

    # determine xyz earth coordinates for lat/long of observatory
    flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
    EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)

    # set earth velocity vectors
    EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])

    # transform into ICRF frame
    EO_bary = sxform("ITRF93", "J2000", epoch) * EO_earth
    @cusync EO_bary_gpu = CuArray(EO_bary)

    # get vector from barycenter to observatory on Earth's surface
    BO_bary = BE_bary .+ EO_bary

    # set string for ltt and abberation
    lt_flag = "CN+S"

    # get light travel time corrected OS vector
    OS_bary, OS_lt, OS_dlt = spkltc(10, epoch, "J2000", lt_flag, BO_bary)
    @cusync OS_bary_gpu = CuArray(OS_bary)

    # get modified epch
    epoch_lt = epoch - OS_lt

    # get rotation matrix for sun
    sun_rot_mat = CuArray(pxform("IAU_SUN", "J2000", epoch_lt))

    threads1 = (16,16)
    blocks1 = cld(prod(CUDA.size(θe)), prod(threads1))
    GRASS.@cusync @cuda threads=threads1 blocks=blocks1 calc_grid_edge_xyz!(ϕe, θe, xs, ys, zs, sun_radius, sun_rot_mat, OS_bary_gpu)
    return Array(xs), Array(ys), Array(zs)
end

function calc_grid_edge_xyz!(ϕe, θe, xs, ys, zs, sun_radius, sun_rot_mat, OS_bary)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = gridDim().y * blockDim().y

    # total number of elements output array
    Nθ_max = CUDA.size(θe, 2)
    num_tiles = prod(CUDA.size(θe))

    for i in idx:sdx:CUDA.length(ϕe)
        for j in idy:sdy:CUDA.size(θe, 2)
            ϕ1 = ϕe[i]
            θ1 = θe[i,j]

            if ((j > 1) & iszero(θ1))
                continue
            end
            
            # get cartesian coords of patch
            x, y, z = sphere_to_cart_gpu_eclipse(sun_radius, ϕ1, θ1)

            # xyz rotated for SP bary
            x_new = sun_rot_mat[1] * x + sun_rot_mat[4] * y + sun_rot_mat[7] * z
            y_new = sun_rot_mat[2] * x + sun_rot_mat[5] * y + sun_rot_mat[8] * z
            z_new = sun_rot_mat[3] * x + sun_rot_mat[6] * y + sun_rot_mat[9] * z

            x_new += OS_bary[1]
            y_new += OS_bary[2]
            z_new += OS_bary[3]

            # set the elements and return
            @inbounds xs[i,j] = x_new
            @inbounds ys[i,j] = y_new
            @inbounds zs[i,j] = z_new
        end 
    end
    return nothing
end

#---------------- here down needs clearing
function calc_europa_quantities_gpu!(epoch::T1, obs_long::T1, obs_lat::T1, alt::T1, wavelength, 
                                      LD_law::String, ext_toggle::T1, ext_coeff::T1,
                                      disk::DiskParamsEclipse{T2}, gpu_allocs::GPUAllocsEclipse{T1}) where {T1<:AF, T2<:AF}

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    projected_v = gpu_allocs.projected_v
    earth_v = gpu_allocs.earth_v
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    ax_codes = gpu_allocs.ax_codes

    # re-zero everything
    @cusync begin
        μs .= 0.0
        ld .= 0.0
        projected_v .= 0.0
        earth_v .= 0.0
        ext .= 0.0
        dA .= 0.0
        z_rot .= 0.0
    end

    # convert scalars from disk params to desired precision
    A = convert(T2, disk.A)
    B = convert(T2, disk.B)
    C = convert(T2, disk.C)
    ext_coeff_gpu = convert(T2, ext_coeff)

    # #europa IR
    # u1 = convert(T2, 0.2692848862319729)
    # u2 = convert(T2, 0.20581805683861285)

    if LD_law == "SSD"
        filtered_df = quad_ld_coeff_SSD[quad_ld_coeff_SSD.wavelength .== wavelength, :]
        u1 = convert(T2, filtered_df.u1[1])
        u2 = convert(T2, filtered_df.u2[1])
    end

    # geometry from disk
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    Nsubgrid = disk.Nsubgrid

    # query SPICE for E, S, M position (km) and velocities (km/s)
    BE_bary = spkssb(399,epoch,"J2000")

    # determine xyz earth coordinates for lat/long of observatory
    flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
    EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)
    # set earth velocity vectors
    EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
    # transform into ICRF frame
    EO_bary = sxform("IAU_EARTH", "J2000", epoch) * EO_earth
    @cusync EO_bary_gpu = CuArray(EO_bary)

    # get vector from barycenter to observatory on Earth's surface
    BO_bary = BE_bary .+ EO_bary

    # set string for ltt and abberation
    lt_flag = "CN+S"

    #Europa position vectors:
    EuE_bary = spkezp(399, epoch, "J2000", lt_flag, 502)[1]
    EuS_bary = spkezp(10, epoch, "J2000", lt_flag, 502)[1]
    phase_angle = acos(calc_mu(EuE_bary, EuS_bary))
    @cusync EuE_bary_gpu = CuArray(EuE_bary)
    @cusync EuS_bary_gpu = CuArray(EuS_bary)
    @cusync moon_radius_gpu = CuArray([moon_radius])
    @cusync earth_radius_gpu = CuArray([earth_radius])
    @cusync sun_radius_gpu = CuArray([sun_radius])

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

    @cusync begin
        Nθ = CuArray{Float64}(disk.Nθ)
    end

    @cusync wavelength_gpu = CuArray(wavelength)
    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(prod(CUDA.size(μs)), prod(threads1))
    @cusync @cuda threads=threads1 blocks=blocks1 calc_europa_quantities_gpu!(wavelength_gpu, μs, z_rot, ax_codes,
                                                                               Nϕ, Nθ, Nsubgrid, Nθ_max, ld, projected_v, earth_v, ext, 
                                                                               dA, moon_radius_gpu, OS_bary_gpu, OM_bary_gpu, EO_bary_gpu,
                                                                               sun_rot_mat_gpu, sun_radius_gpu, A, B, C, u1, u2, 
                                                                               ext_toggle, ext_coeff_gpu, earth_radius_gpu, CB1, CB2, CB3, EuE_bary_gpu, EuS_bary_gpu, CuArray(collect(theta_values)), CuArray(SH_y), CuArray(CB_y))
    return
end                               

function calc_europa_quantities_gpu!(wavelength, μs, z_rot, ax_codes,
                                      Nϕ, Nθ, Nsubgrid, Nθ_max, ld, projected_v, earth_v, ext,
                                      dA, moon_radius, OS_bary, OM_bary, EO_bary,
                                      sun_rot_mat, sun_radius, A, B, C, u1, u2, 
                                      ext_toggle, ext_coeff_gpu, earth_radius, CB1, CB2, CB3, EuE_bary_gpu, EuS_bary_gpu, theta_values, SH_y, CB_y) 
    sun_radius = sun_radius[1] 
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
        z_rot_numerator = CUDA.zero(CUDA.eltype(μs))
        z_rot_denominator = CUDA.zero(CUDA.eltype(μs))
        dA_sum = CUDA.zero(CUDA.eltype(μs))
        ld_sum = CUDA.zero(CUDA.eltype(μs))
        projected_v_sum = CUDA.zero(CUDA.eltype(μs))
        earth_v_sum = CUDA.zero(CUDA.eltype(μs))
        ext_sum = CUDA.zero(CUDA.eltype(μs))
        x_sum = CUDA.zero(CUDA.eltype(μs))
        y_sum = CUDA.zero(CUDA.eltype(μs))
        z_sum = CUDA.zero(CUDA.eltype(μs))

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

                # get cartesian coords of patch
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

                # xyz rotated for SP bary
                x_new = sun_rot_mat[1] * x + sun_rot_mat[4] * y + sun_rot_mat[7] * z
                y_new = sun_rot_mat[2] * x + sun_rot_mat[5] * y + sun_rot_mat[8] * z
                z_new = sun_rot_mat[3] * x + sun_rot_mat[6] * y + sun_rot_mat[9] * z

                # vel component rotated for SP bary
                vx = sun_rot_mat[1] * d + sun_rot_mat[4] * e + sun_rot_mat[7] * f
                vy = sun_rot_mat[2] * d + sun_rot_mat[5] * e + sun_rot_mat[8] * f
                vz = sun_rot_mat[3] * d + sun_rot_mat[6] * e + sun_rot_mat[9] * f

                # OP_bary state vector
                OP_bary_x = OS_bary[1] + x_new
                OP_bary_y = OS_bary[2] + y_new
                OP_bary_z = OS_bary[3] + z_new

                # Europa to position state vector
                EuP_bary_x = EuS_bary_gpu[1] + x_new
                EuP_bary_y = EuS_bary_gpu[2] + y_new
                EuP_bary_z = EuS_bary_gpu[3] + z_new

                # calculate mu
                μ_sub = calc_mu_gpu(x_new, y_new, z_new, OP_bary_x, OP_bary_y, OP_bary_z) 
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub
                μ_count += 1

                # sum on vector components
                x_sum += x
                y_sum += y
                z_sum += z

                # get OP_bary and SP_bary between them and find projected_velocities_no_cb
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n2 = CUDA.sqrt(vx^2.0 + vy^2.0 + vz^2.0)
                angle = (OP_bary_x * vx + OP_bary_y * vy + OP_bary_z * vz) / (n1 * n2)
                # v_rot_sub = (n2 * angle) + (CB1*CUDA.exp(μ_sub)^2.0 + CB2*CUDA.exp(μ_sub) + CB3)
                v_rot_sub = (n2 * angle)
                v_rot_sub *= 1000.0

                n2 = CUDA.sqrt(OS_bary[4]^2.0 + OS_bary[5]^2.0 + OS_bary[6]^2.0)
                angle = (OP_bary_x * OS_bary[4] + OP_bary_y * OS_bary[5] + OP_bary_z * OS_bary[6]) / (n1 * n2)
                v_orbit_sub = (n2 * angle)
                v_orbit_sub *= 1000.0

                # get projected area element
                dA_sub = calc_dA_gpu(sun_radius, ϕc_sub, dϕ, dθ)
                dA_sub *= μ_sub
                dA_sum += dA_sub

                # calculate distance (Europa)
                n1 = CUDA.sqrt(EuP_bary_x^2.0 + EuP_bary_y^2.0 + EuP_bary_z^2.0)  
                n2 = CUDA.sqrt(EuE_bary_gpu[1]^2.0 + EuE_bary_gpu[2]^2.0 + EuE_bary_gpu[3]^2.0)  
                d2 = acos((EuE_bary_gpu[1] * EuP_bary_x + EuE_bary_gpu[2] * EuP_bary_y + EuE_bary_gpu[3] * EuP_bary_z) / (n2 * n1))
               
                z_rot_sub = v_rot_sub / c_ms
                z_rot_sub += (v_orbit_sub / c_ms)

                projected_v_sum += v_rot_sub
                earth_v_sum += v_orbit_sub

                # iterate counter 
                count += 1

                # zenith
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n3 = CUDA.sqrt(EO_bary[1]^2.0 + EO_bary[2]^2.0 + EO_bary[3]^2.0) 
                zenith = (CUDA.acos((EO_bary[1] * OP_bary_x + EO_bary[2] * OP_bary_y + EO_bary[3] * OP_bary_z) / (n3 * n1)))

                for wl in eachindex(wavelength)
                    # get limb darkening
                    ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2)

                    if (rad2deg(d2) < 0.4)
                        # continue
                        ld_sum += ld_sub*lorentzian_phase_curve(d2)
                    else 
                        ld_sum += ld_sub
                    end

                    if CUDA.isone(ext_toggle)
                        ext_sub = CUDA.exp(-((1/(CUDA.cos(zenith)))*ext_coeff_gpu))
                        ext_sum += ext_sub
                    end
                end

                # recalculate value since local vars from wl for loop are lost
                ld_sub = quad_limb_darkening_gpu(μ_sub, u1, u2)
                ext_sub = CUDA.exp(-((1/(CUDA.cos(zenith)))*ext_coeff_gpu))

                if ext_toggle == 1.0
                    z_rot_numerator += z_rot_sub * dA_sub * ld_sub * ext_sub
                    z_rot_denominator += dA_sub * ld_sub * ext_sub
                else
                    z_rot_numerator += z_rot_sub * dA_sub * ld_sub
                    z_rot_denominator += dA_sub * ld_sub 
                end
            end
        end
        # take averages
        @inbounds μs[m,n] = μ_sum / μ_count
        @inbounds dA[m,n] = dA_sum 

        if count > 0
            # set scalar quantity elements as average
            @inbounds projected_v[m,n] = projected_v_sum / count
            @inbounds earth_v[m,n] = earth_v_sum / count

            @inbounds z_rot[m,n] = z_rot_numerator / z_rot_denominator

            for wl in eachindex(wavelength)
                @inbounds ld[m,n,wl] = ld_sum / count
                @inbounds ext[m,n,wl] = ext_sum / count
            end

            # set vector components as average
            @inbounds xx = x_sum / μ_count
            @inbounds yy = y_sum / μ_count
            @inbounds zz = z_sum / μ_count

            @inbounds ax_codes[m, n] = find_nearest_ax_gpu(yy / sun_radius, zz / sun_radius)
        end
    end

    return nothing
end