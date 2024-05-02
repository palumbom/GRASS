function calc_eclipse_quantities_gpu!(t::Int, epoch, obs_long, obs_lat, alt, wavelength,
                                        disk::DiskParams{T2}, gpu_allocs::GPUAllocs{T1},
                                        eclipse_allocs::EclipseAllocs{T2}) where {T1<:AF, T2<:AF}

    lwavgrid = eclipse_allocs.lwavgrid
    rwavgrid = eclipse_allocs.rwavgrid
    allwavs = eclipse_allocs.allwavs
    allints = eclipse_allocs.allints
    bist = eclipse_allocs.bist
    intt = eclipse_allocs.intt
    widt = eclipse_allocs.widt
    ϕc = eclipse_allocs.ϕc
    θc = eclipse_allocs.θc
    μs = eclipse_allocs.μs
    ld = eclipse_allocs.ld
    ext = eclipse_allocs.ext
    dA = eclipse_allocs.dA
    wts = eclipse_allocs.wts
    xyz = eclipse_allocs.xyz
    cbs = eclipse_allocs.cbs
    z_rot = eclipse_allocs.z_rot
    ax_codes = eclipse_allocs.ax_codes
    keys = eclipse_allocs.keys
    dA_total_proj_mean = eclipse_allocs.dA_total_proj_mean
    mean_intensity = eclipse_allocs.mean_intensity
    mean_weight_v_no_cb = eclipse_allocs.mean_weight_v_no_cb
    mean_weight_v_earth_orb = eclipse_allocs.mean_weight_v_earth_orb
    pole_vector_grid = eclipse_allocs.pole_vector_grid
    SP_sun_pos = eclipse_allocs.SP_sun_pos
    SP_sun_vel = eclipse_allocs.SP_sun_vel
    SP_bary = eclipse_allocs.SP_bary
    SP_bary_pos = eclipse_allocs.SP_bary_pos
    SP_bary_vel = eclipse_allocs.SP_bary_vel
    OP_bary = eclipse_allocs.OP_bary
    mu_grid = eclipse_allocs.mu_grid
    projected_velocities_no_cb = eclipse_allocs.projected_velocities_no_cb
    distance = eclipse_allocs.distance
    v_scalar_grid = eclipse_allocs.v_scalar_grid
    v_earth_orb_proj = eclipse_allocs.v_earth_orb_proj

    # convert scalars from disk params to desired precision
    sun_radius = convert(T2, disk.sun_radius)
    A = convert(T2, disk.A)
    B = convert(T2, disk.B)
    C = convert(T2, disk.C)
    v0 = convert(T2, disk.v0)
    u1 = convert(T2, disk.u1)
    u2 = convert(T2, disk.u2)

    # geometry from disk
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    Nsubgrid = disk.Nsubgrid
    Nϕ_sub = Nϕ * Nsubgrid
    Nθ_sub = maximum(disk.Nθ) * Nsubgrid

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

    # get vector from observatory on earth's surface to moon center
    OM_bary, OM_lt, OM_dlt = spkltc(301, epoch, "J2000", lt_flag, BO_bary)

    # get modified epch
    epoch_lt = epoch - OS_lt

    # get rotation matrix for sun
    sun_rot_mat = pxform("IAU_SUN", "J2000", epoch_lt)

    # compute geometric parameters, average over subtiles
    threads1 = 256
    blocks1 = cld(CUDA.length(μs), prod(threads1))
    @cusync @captured @cuda threads=threads1 blocks=blocks1 calc_eclipse_quantities_gpu!(t, wavelength, ϕc, θc, μs, wts, z_rot,
                                                                                          Nϕ, Nsubgrid, ld, ext, dA, xyz, ax_codes,
                                                                                          dA_total_prof_mean, mean_intensity, mean_weight_v_no_cb, 
                                                                                          mean_weight_v_earth_orb, pole_vector_grid,
                                                                                          SP_sun_pos, SP_sun_vel, SP_bary, SP_bary_pos, 
                                                                                          SP_bary_vel, OP_bary, mu_grid, projected_velocities_no_cb,
                                                                                          distance, v_scalar_grid, v_earth_orb_proj,
                                                                                          moon_radius, EO_bary, OS_bary, OM_bary, sun_rot_mat, 
                                                                                          sun_radius, 
                                                                                          R_x, O⃗, A, B, C, v0, u1, u2)
    return nothing
end                               

function calc_eclipse_quantities_gpu!(t, wavelength, ϕc, θc, μs, wts, z_rot,
                                        Nϕ, Nsubgrid, ld, ext, dA, xyz, ax_codes,
                                        dA_total_prof_mean, mean_intensity, mean_weight_v_no_cb, 
                                        mean_weight_v_earth_orb, pole_vector_grid,
                                        SP_sun_pos, SP_sun_vel, SP_bary, SP_bary_pos, 
                                        SP_bary_vel, OP_bary, mu_grid, projected_velocities_no_cb,
                                        distance, v_scalar_grid, v_earth_orb_proj,
                                        moon_radius, EO_bary, OS_bary, OM_bary, sun_rot_mat, 
                                        sun_radius, 
                                        R_x, O⃗, A, B, C, v0, u1, u2)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    # get latitude subtile step size
    N_ϕ_edges = Nϕ * Nsubgrid
    dϕ = (CUDA.deg2rad(90.0) - CUDA.deg2rad(-90.0)) / (N_ϕ_edges)

    # linear index over course grid tiles
    for i in idx:sdx:CUDA.length(μs)
        # get dϕ and dθ for large tiles
        dϕ_large = (CUDA.deg2rad(90.0) - CUDA.deg2rad(-90.0)) / Nϕ
        dθ_large = (CUDA.deg2rad(360.0) - CUDA.deg2rad(0.0)) / get_Nθ(ϕc[i], dϕ_large)

        # get starting edge of large tiles
        ϕl = ϕc[i] - dϕ_large / 2.0
        θl = θc[i] - dθ_large / 2.0

        # get number of longitude tiles in course latitude slice
        N_θ_edges = get_Nθ(ϕc[i], dϕ_large)
        dθ = (CUDA.deg2rad(360.0) - CUDA.deg2rad(0.0)) / (N_θ_edges * Nsubgrid)

        # set up sum holders for scalars
        μ_sum = CUDA.zero(CUDA.eltype(μs))
        v_sum = CUDA.zero(CUDA.eltype(z_rot))
        v_sum_earth = CUDA.zero(CUDA.eltype(z_rot))
        ld_sum = CUDA.zero(CUDA.eltype(wts))
        dA_sum = CUDA.zero(CUDA.eltype(wts))
        ext_sum = CUDA.zero(CUDA.eltype(wts))

        # initiate counter
        μ_count = 0
        count = 0

        # loop over ϕ subgrid
        for j in 1:Nsubgrid
            # get coordinates of latitude subtile center
            ϕsub = ϕl + (dϕ/2.0) + (j - 1) * dϕ

            # loop over θsubgrid
            for k in 1:Nsubgrid
                # get coordinates of longitude subtile center
                θsub = θl + (dθ/2.0) + (k - 1) * dθ

                # get cartesian coords of patch
                x, y, z = sphere_to_cart_gpu_eclipse(sun_radius, ϕsub, θsub)

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
                rp = 2π * sun_radius * CUDA.cos(ϕsub) / rotation_period_gpu(ϕsub, A, B, C)

                # get in units of c
                rp /= 86400.0

                # set magnitude of vector
                d *= rp
                e *= rp
                f *= rp

                #xyz rotated for SP bary
                x = sun_rot_mat[1] * a + sun_rot_mat[2] * a + sun_rot_mat[3] * a
                y = sun_rot_mat[4] * b + sun_rot_mat[5] * b + sun_rot_mat[6] * b
                z = sun_rot_mat[7] * c + sun_rot_mat[8] * c + sun_rot_mat[9] * c
                #vel component rotated for SP bary
                vx = sun_rot_mat[1] * d + sun_rot_mat[2] * d + sun_rot_mat[3] * d
                vy = sun_rot_mat[4] * e + sun_rot_mat[5] * e + sun_rot_mat[6] * e
                vz = sun_rot_mat[7] * f + sun_rot_mat[8] * f + sun_rot_mat[9] * f

                #OP_bary state vector
                OP_bary_x = OS_bary[1] + x
                OP_bary_y = OS_bary[2] + y
                OP_bary_z = OS_bary[3] + z

                OP_bary_vx = OS_bary[4] + vx
                OP_bary_vy = OS_bary[5] + vy
                OP_bary_vz = OS_bary[6] + vz

                # calculate mu
                μ_sub = calc_mu_gpu(x, y, z, O⃗)
                if μ_sub <= 0.0
                    continue
                end
                μ_sum += μ_sub
                μ_count += 1

                # get OP_bary and SP_bary between them and find projected_velocities_no_cb
                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)
                n2 = CUDA.sqrt(vx^2.0 + vy^2.0 + vz^2.0)
                angle = (OP_bary_x * vx + OP_bary_y * vy + OP_bary_z * vz) / (n1 * n2)
                v_sum += (n2 * angle)
                v_sum *= 1000.0
                z_rot_sub = v_sum / c_ms

                n1 = CUDA.sqrt(OP_bary_x^2.0 + OP_bary_y^2.0 + OP_bary_z^2.0)  
                n2 = CUDA.sqrt(OS_bary[4]^2.0 + OS_bary[5]^2.0 + OS_bary[6]^2.0)
                angle = (OP_bary_x * OS_bary[4] + OP_bary_y * OS_bary[5] + OP_bary_z * OS_bary[6]) / (n1 * n2)
                v_sum_earth += (n2 * angle)
                v_sum_earth *= 1000.0
                z_rot_sub += (v_sum / c_ms)

                #calculate distance
                n1 = CUDA.sqrt(OM_bary[1]^2.0 + OM_bary[2]^2.0 + OM_bary[3]^2.0)  
                n2 = CUDA.sqrt(OP_bary_x[4]^2.0 + OP_bary_y[5]^2.0 + OP_bary_z[6]^2.0)
                d2 = (OM_bary[1] * OP_bary_x + OM_bary[2] * OP_bary_y + OM_bary[3] * OP_bary_z) / (n1 * n2)
                if (d2 <= rad_planet^2.0) & (z_planet > 0.0)
                    continue
                end

#---------------------

                # get limb darkening
                ld = quad_limb_darkening_gpu(μ_sub, u1, u2)
                ld_sum += ld

                # get projected area element
                dA = calc_dA_gpu(sun_radius, ϕsub, dϕ, dθ)
                dA *= μ_sub
                dA_sum += dA

                # iterate counter
                count += 1
            end
        end

        # take averages
        if count > 0
            # set scalar quantity elements as average
            @inbounds μs[i] = μ_sum / μ_count
            @inbounds z_rot[i] = v_sum / ld_sum
            @inbounds wts[i] = (ld_sum / count) * dA_sum
        else
            @inbounds μs[i] = 0.0
            @inbounds z_rot[i] = 0.0
            @inbounds wts[i] = 0.0
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