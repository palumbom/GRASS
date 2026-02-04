function eclipse_compute_quantities!(epoch::T, t, obs_long::T, obs_lat::T, alt::T, wavelength, 
                                      LD_law::String, ext_toggle, ext_coeff::T,
                                      disk::DiskParamsEclipse{T}, cpu_allocs::SynthWorkspaceEclipse{T}) where {T<:AF}

    μs = cpu_allocs.μs
    ld = cpu_allocs.ld
    ext = cpu_allocs.ext
    dA = cpu_allocs.dA
    z_rot = cpu_allocs.z_rot
    ax_codes = cpu_allocs.ax_codes
    SP_sun_pos = cpu_allocs.SP_sun_pos
    SP_sun_vel = cpu_allocs.SP_sun_vel
    v_scalar_grid = cpu_allocs.v_scalar_grid
    pole_vector_grid = cpu_allocs.pole_vector_grid
    SP_bary_pos = cpu_allocs.SP_bary_pos
    SP_bary_vel = cpu_allocs.SP_bary_vel
    mu_grid = cpu_allocs.mu_grid
    SP_bary = cpu_allocs.SP_bary
    OP_bary = cpu_allocs.OP_bary
    projected_velocities_no_cb = cpu_allocs.projected_velocities_no_cb
    v_earth_orb_proj = cpu_allocs.v_earth_orb_proj
    distance = cpu_allocs.distance
    xyz = cpu_allocs.xyz
    mean_weight_v_no_cb = cpu_allocs.mean_weight_v_no_cb
    mean_weight_v_earth_orb = cpu_allocs.mean_weight_v_earth_orb

    A = disk.A
    B = disk.B
    C = disk.C

    if LD_law == "300_quadratic"
        filtered_df = quad_ld_coeff_300[quad_ld_coeff_300.wavelength .== wavelength, :]
        u1 = filtered_df.u1[1]
        u2 = filtered_df.u2[1]
        u3 = NaN
        u4 = NaN
    end

    if LD_law == "300_4parameter"
        filtered_df = quad_ld_coeff_300[quad_ld_coeff_300.wavelength .== wavelength, :]
        u1 = filtered_df.law_4u1[1]
        u2 = filtered_df.law_4u2[1]
        u3 = filtered_df.law_4u3[1]
        u4 = filtered_df.law_4u4[1]
    end

    if LD_law == "HD_quadratic"
        filtered_df = quad_ld_coeff_HD[quad_ld_coeff_HD.wavelength .== wavelength, :]
        u1 = filtered_df.u1[1]
        u2 = filtered_df.u2[1]
        u3 = NaN
        u4 = NaN
    end

    if LD_law == "HD_4parameter"
        filtered_df = quad_ld_coeff_HD[quad_ld_coeff_HD.wavelength .== wavelength, :]
        u1 = filtered_df.law_4u1[1]
        u2 = filtered_df.law_4u2[1]
        u3 = filtered_df.law_4u3[1]
        u4 = filtered_df.law_4u4[1]
    end

    if LD_law == "SSD_quadratic"
        filtered_df = quad_ld_coeff_SSD[quad_ld_coeff_SSD.wavelength .== wavelength, :]
        u1 = filtered_df.u1[1]
        u2 = filtered_df.u2[1]
        u3 = NaN
        u4 = NaN
    end

    if LD_law == "SSD_4parameter"
        filtered_df = quad_ld_coeff_SSD[quad_ld_coeff_SSD.wavelength .== wavelength, :]
        u1 = filtered_df.law_4u1[1]
        u2 = filtered_df.law_4u2[1]
        u3 = filtered_df.law_4u3[1]
        u4 = filtered_df.law_4u4[1]
    end

    # geometry from disk
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    Nsubgrid = disk.Nsubgrid

    # query for E, S, M position (km) and velocities (km/s)
    BE_bary = spkssb(399,epoch,"J2000")

    # determine xyz earth coordinates for lat/long of observatory
    flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
    EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)
    #set earth velocity vectors
    EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
    # transform into ICRF frame
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

    # loop over disk positions
    for i in eachindex(disk.ϕc)
        for j in 1:disk.Nθ[i]
            # subdivide the tile
            ϕe_sub = range(disk.ϕe[i], disk.ϕe[i+1], length=Nsubgrid+1)
            θe_sub = range(disk.θe[i,j], disk.θe[i,j+1], length=Nsubgrid+1)
            ϕc_sub = get_grid_centers(ϕe_sub)
            θc_sub = get_grid_centers(θe_sub)
            subgrid = Iterators.product(ϕc_sub, θc_sub)

            # get cartesian coord for each subgrid 
            SP_sun_pos .= map(x -> pgrrec("SUN", getindex(x,2), getindex(x,1), 0.0, sun_radius, 0.0), subgrid)

            # get differential rotation velocities
            v_scalar_grid .= map(x -> v_scalar(x...), subgrid)
            # convert v_scalar to from km/day km/s
            v_scalar_grid ./= 86400.0

            # determine pole vector for each patch
            pole_vector_grid!(SP_sun_pos, pole_vector_grid)

            # get velocity vector direction and set magnitude
            v_vector(SP_sun_pos, pole_vector_grid, v_scalar_grid, SP_sun_vel)

            for k in eachindex(SP_sun_pos)
                SP_bary_pos[k] .= (sun_rot_mat * SP_sun_pos[k])
                SP_bary_vel[k] .= (sun_rot_mat * SP_sun_vel[k])
                SP_bary[k] = vcat(SP_bary_pos[k], SP_bary_vel[k])
            end

            # get vector from obs to each patch on Sun's surface
            for k in eachindex(OP_bary)
                OP_bary[k] = OS_bary .+ SP_bary[k]
            end

            # calculate mu at each point
            calc_mu_grid!(SP_bary, OP_bary, mu_grid)
            # move on if everything is off the grid
            all(mu_grid .<= zero(T)) && continue

            # get projected velocity for each patch
            projected!(SP_bary, OP_bary, projected_velocities_no_cb)
            # convert from km/s to m/s
            projected_velocities_no_cb .*= 1000.0
            z_rot_sub = (projected_velocities_no_cb)./c_ms

            # get relative orbital motion in m/s
            v_delta = OS_bary[4:6]
            for k in eachindex(v_earth_orb_proj)
                B = view(OP_bary[k], 1:3)
                angle = dot(B, v_delta) / (norm(B) * norm(v_delta))
                v_earth_orb_proj[k] = norm(v_delta) * angle
            end
            # convert from km/s to m/s
            v_earth_orb_proj .*= 1000.0
            z_rot_sub .+= v_earth_orb_proj/c_ms

            # calculate the distance between tile corner and moon
            for i in eachindex(OP_bary)
                distance[i] = calc_proj_dist(OM_bary[1:3], OP_bary[i][1:3])
            end

            # get indices for visible patches
            idx1 = mu_grid .> 0.0
            idx3 = (idx1) .& (distance .> atan((moon_radius)/norm(OM_bary[1:3]))) 
            # assign the mean mu as the mean of visible mus
            μs[i,j,t] = mean(view(mu_grid, idx1))

            # find xz at mean value of mu and get axis code (i.e., N, E, S, W)
            xyz[i,j,1] = mean(view(getindex.(SP_sun_pos,1), idx1))
            xyz[i,j,2] = mean(view(getindex.(SP_sun_pos,2), idx1))
            xyz[i,j,3] = mean(view(getindex.(SP_sun_pos,3), idx1)) 
            if xyz[i,j,2] !== NaN
                ax_codes[i,j,t] = find_nearest_ax_code_eclipse(xyz[i,j,2]/sun_radius, xyz[i,j,3]/sun_radius) 
            end
            
            # calculate area element of tile
            dϕ = step(ϕe_sub) 
            dθ = step(θe_sub) 
            dA_sub = map(x -> calc_dA(1.0, getindex(x,1), dϕ, dθ), subgrid)

            # get total projected, visible area of larger tile
            dA_total_proj = dA_sub .* mu_grid
            dA[i,j,t] = sum(view(dA_total_proj, idx1))   

            mean_weight_v_no_cb[i,j,t] = mean(view(projected_velocities_no_cb, idx3))
            mean_weight_v_earth_orb[i,j,t] = mean(view(v_earth_orb_proj, idx3))

            for l in eachindex(wavelength)
                # limb darkening
                if isnan(u3)
                    ld_sub = map(x -> quad_limb_darkening(x, u1, u2), mu_grid)
                else
                    ld_sub = map(x -> quad_limb_darkening(x, u1, u2, u3, u4), mu_grid)
                end

                ld[i,j,l] = mean(view(ld_sub, idx3))  

                # z_rot
                if ext_toggle == false
                    z_rot[i,j,l] = sum(view(z_rot_sub .* ld_sub .* dA_total_proj, idx3)) ./ sum(view(ld_sub .* dA_total_proj, idx3))
                end

                # z_rot + extinction
                if ext_toggle == true
                    extin = map(x -> exp(-((1/cos(x))*ext_coeff)), map(x -> calc_proj_dist(x[1:3], EO_bary[1:3]), OP_bary))
                    ext[i,j,l] = mean(view(extin, idx3)) 

                    z_rot[i,j,l] = sum(view(z_rot_sub .* ld_sub .* dA_total_proj .* extin, idx3)) ./ sum(view(ld_sub .* dA_total_proj .* extin, idx3))

                    if isnan(ext[i,j,l])
                        ld[i,j,l] = 0.0
                        ext[i,j,l] = 0.0
                    end                
                end
       
                if isnan(ld[i,j,l])
                    ld[i,j,l] = 0.0
                    ext[i,j,l] = 0.0
                    mean_weight_v_no_cb[i,j,l] = 0.0
                    mean_weight_v_earth_orb[i,j,l] = 0.0
                end
            end
        end
    end
    return 
end

function generate_tloop_eclipse!(tloop::AA{Int}, wsp::SynthWorkspaceEclipse{T}, soldata::SolarData{T}, t) where T<:AF
    generate_tloop_eclipse!(tloop, wsp.μs[:, :, t], wsp.keys, soldata.len)
    return nothing
end

function generate_tloop_eclipse!(tloop::AA{Int}, μs::AA{T}, keys::AA{Tuple{Symbol, Symbol}},
                         len::Dict{Tuple{Symbol, Symbol}, Int64}) where T<:AF
    # loop over μ positions
    for i in eachindex(μs)
        # move on if we are off the grid
        μs[i] <= zero(T) && continue

        # get length of input data for place on disk
        maxlen = len[keys[i]]

        # generate random index
        tloop[i] = floor(Int, rand() * maxlen) + 1
    end
    return nothing
end

function get_keys_and_cbs_eclispe!(wsp::SynthWorkspaceEclipse{T}, soldata::SolarData{T}, t) where T<:AF
    get_keys_and_cbs_eclispe!(wsp.keys, wsp.μs[:, :, t], wsp.cbs, wsp.ax_codes[:, :, t], soldata)
    return nothing
end


function get_keys_and_cbs_eclispe!(keys::AA{Tuple{Symbol, Symbol}}, μs::AA{T}, cbs::AA{T},
                           ax_codes::AA{Int}, soldata::SolarData{T2}) where  {T<:AF, T2<:AF}
    # get the mu and axis codes
    disc_mu = soldata.mu
    disc_ax = soldata.ax

    # type finagling
    disc_mu = convert.(T, disc_mu)

    # loop over μ positions
    for i in eachindex(μs)
        # move on if we are off the grid
        μs[i] <= zero(T) && continue

        # get input data for place on disk
        the_key = get_key_for_pos(μs[i], ax_codes[i], disc_mu, disc_ax)
        cbs[i] = convert(T, soldata.cbs[the_key])
        keys[i] = the_key
    end
    return nothing
end

function calc_rossiter_quantities!(xyz_planet::AA{T,1}, planet::Planet{T},
                                   disk::DiskParams{T}, wsp::SynthWorkspace{T},
                                   ros_allocs::RossiterAllocs{T}) where T<:AF
    # parse out composite type fields
    Nsubgrid = disk.Nsubgrid

    # allocate memory that wont be needed outside this function
    d2_sub = ros_allocs.d2_sub
    μs_sub = ros_allocs.μs_sub
    ld_sub = ros_allocs.ld_sub
    dA_sub = ros_allocs.dA_sub
    dp_sub = ros_allocs.dp_sub
    xyz_sub = repeat([zeros(3)], Nsubgrid, Nsubgrid)
    z_rot_sub = ros_allocs.z_rot_sub
    idx1 = ros_allocs.idx1
    idx2 = ros_allocs.idx2
    idx3 = ros_allocs.idx3

    # loop over disk positions
    for i in eachindex(wsp.ϕc)
        # get the number of theta tiles needed for the latitude tiles
        Nθ = get_Nθ(wsp.ϕc[i], step(disk.ϕe))

        # get the latitude and longitude increments
        dϕ = step(disk.ϕe)
        dθ = deg2rad(360.0) / Nθ

        # get edges of large tile
        ϕ_l = wsp.ϕc[i] - dϕ/2.0
        ϕ_r = wsp.ϕc[i] + dϕ/2.0
        θ_l = wsp.θc[i] - dθ/2.0
        θ_r = wsp.θc[i] + dθ/2.0

        # subdivide the tile
        ϕe_sub = range(ϕ_l, ϕ_r, length=Nsubgrid+1)
        θe_sub = range(θ_l, θ_r, length=Nsubgrid+1)
        ϕc_sub = get_grid_centers(ϕe_sub)
        θc_sub = get_grid_centers(θe_sub)
        subgrid = Iterators.product(ϕc_sub, θc_sub)

        # get cartesian coord for each subgrid and rotate by rot. matrix
        xyz_sub .= map(x -> sphere_to_cart.(disk.ρs, x...), subgrid)
        xyz_sub .= map(x -> disk.R_x * x, xyz_sub)

        # calculate mu at each point
        μs_sub .= map(x -> calc_mu(x, disk.O⃗), xyz_sub)
        if all(μs_sub .<= zero(T))
            continue
        end

        # calculate the distance between subtile center and planet
        d2_sub .= map(x -> calc_proj_dist2(x, xyz_planet), xyz_sub)

        # if entire course tile visible, use old weights and move on
        if all(d2_sub .> planet.radius^2.0) | (xyz_planet[3] < 0.0)
            continue
        end

        # calc limb darkening
        ld_sub .= map(x -> quad_limb_darkening(x, disk.u1, disk.u2), μs_sub)

        # get rotational velocity for location on disk
        z_rot_sub .= map(x -> patch_velocity_los(x..., disk), subgrid)

        # calculate area element of tile
        dA_sub .= map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕe_sub), step(θe_sub)), subgrid)
        dp_sub .= map(x -> abs(dot(x .- disk.O⃗, x)), xyz_sub) / norm(disk.O⃗)

        # get indices for visible patches
        idx1 .= μs_sub .> 0.0
        idx2 .= d2_sub .> planet.radius^2.0
        idx3 .= idx1 .& idx2

        # if no patches are visible, set wts, etc. to zero and move on
        if all(iszero(idx3))
            ros_allocs.μs[i] = 0.0
            ros_allocs.ld[i] = 0.0
            ros_allocs.dA[i] = 0.0
            ros_allocs.wts[i] = 0.0
            ros_allocs.z_rot[i] = 0.0
            continue
        end

        # get total projected, visible area of larger tile
        dA_total = sum(view(dA_sub, idx3))
        dA_total_proj = sum(view(dA_sub .* dp_sub, idx3))

        # set limb darkening as mean of visible patches
        ros_allocs.μs[i] = mean(view(μs_sub, idx1))
        ros_allocs.ld[i] = mean(view(ld_sub, idx3))
        ros_allocs.dA[i] = dA_total_proj

        ros_allocs.wts[i] = mean(view(ld_sub .* dA_total_proj, idx3))
        ros_allocs.z_rot[i] = sum(view(z_rot_sub .* ld_sub, idx3)) ./ sum(view(ld_sub, idx3))
    end
    return nothing
end
