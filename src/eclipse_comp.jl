function eclipse_compute_quantities!(disk::DiskParamsEclipse{T}, epoch::T, obs_long::T,
                                     obs_lat::T, alt::T, wavelength::Vector{Float64}, ϕc::AA{T,2}, θc::AA{T,2},
                                     μs::AA{T,2}, ld::AA{T,3}, ext::AA{T,3} , dA::AA{T,2},
                                     xyz::AA{T,3}, z_rot::AA{T,3}, ax_codes::AA{Int64, 2}, 
                                     dA_total_proj_mean::AA{T,2}, mean_intensity::AA{T,3}, mean_weight_v_no_cb::AA{T,2},
                                     mean_weight_v_earth_orb::AA{T,2}, pole_vector_grid::Matrix{Vector{Float64}},
                                     SP_sun_pos::Matrix{Vector{Float64}}, SP_sun_vel::Matrix{Vector{Float64}}, SP_bary::Matrix{Vector{Float64}}, SP_bary_pos::Matrix{Vector{Float64}},
                                     SP_bary_vel::Matrix{Vector{Float64}}, OP_bary::Matrix{Vector{Float64}}, mu_grid::AA{T,2}, projected_velocities_no_cb::AA{T,2}, 
                                     distance::AA{T,2}, v_scalar_grid::AA{T,2}, v_earth_orb_proj::Matrix{Float64}) where T<:AF

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

    # parse out composite type fields
    Nsubgrid = disk.Nsubgrid

    # loop over disk positions
    for i in eachindex(disk.ϕc)
        for j in 1:disk.Nθ[i]
            # save the tile position
            ϕc[i,j] = disk.ϕc[i]
            θc[i,j] = disk.θc[i,j]

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
            #convert v_scalar to from km/day km/s
            v_scalar_grid ./= 86400.0

            #determine pole vector for each patch
            pole_vector_grid!(SP_sun_pos, pole_vector_grid)

            #get velocity vector direction and set magnitude
            v_vector(SP_sun_pos, pole_vector_grid, v_scalar_grid, SP_sun_vel)

            for k in eachindex(SP_sun_pos)
                SP_bary_pos[k] .= (sun_rot_mat * SP_sun_pos[k])
                SP_bary_vel[k] .= (sun_rot_mat * SP_sun_vel[k])
                SP_bary[k] = vcat(SP_bary_pos[k], SP_bary_vel[k])
            end

            #get vector from obs to each patch on Sun's surface
            for k in eachindex(OP_bary)
                OP_bary[k] = OS_bary .+ SP_bary[k]
            end

            # calculate mu at each point
            calc_mu_grid!(SP_bary, OP_bary, mu_grid)
            # move on if everything is off the grid
            all(mu_grid .< zero(T)) && continue

            #get projected velocity for each patch
            projected!(SP_bary, OP_bary, projected_velocities_no_cb)
            # convert from km/s to m/s
            projected_velocities_no_cb .*= 1000.0
            z_rot_sub = (projected_velocities_no_cb)./c_ms

            #get relative orbital motion in m/s
            v_delta = OS_bary[4:6]
            for k in eachindex(v_earth_orb_proj)
                B = view(OP_bary[k], 1:3)
                angle = dot(B, v_delta) / (norm(B) * norm(v_delta))
                v_earth_orb_proj[k] = norm(v_delta) * angle
            end
            # convert from km/s to m/s
            v_earth_orb_proj .*= 1000.0

            z_rot_sub .+= v_earth_orb_proj/c_ms

            #determine patches that are blocked by moon 
            #calculate the distance between tile corner and moon
            for i in eachindex(OP_bary)
                distance[i] = calc_proj_dist(OM_bary[1:3], OP_bary[i][1:3])
            end

            #get indices for visible patches
            idx1 = mu_grid .> 0.0

            idx3 = (idx1) .& (distance .> atan((moon_radius)/norm(OM_bary[1:3])))

            # assign the mean mu as the mean of visible mus
            μs[i,j] = mean(view(mu_grid, idx1))

            # find xz at mean value of mu and get axis code (i.e., N, E, S, W)
            xyz[i,j,1] = mean(view(getindex.(SP_sun_pos,1), idx1))
            xyz[i,j,2] = mean(view(getindex.(SP_sun_pos,2), idx1))
            xyz[i,j,3] = mean(view(getindex.(SP_sun_pos,3), idx1)) 
            if xyz[i,j,2] !== NaN
                ax_codes[i,j] = find_nearest_ax_code_eclipse(xyz[i,j,2]/sun_radius, xyz[i,j,3]/sun_radius) 
            end
            
            # calculate area element of tile
            dϕ = step(ϕe_sub) 
            dθ = step(θe_sub) 
            dA_sub = map(x -> calc_dA(sun_radius, getindex(x,1), dϕ, dθ), subgrid)
            #get total projected, visible area of larger tile
            dA_total_proj = dA_sub .* mu_grid
            dA_total_proj_mean[i,j] = sum(view(dA_total_proj, idx1))    
                 
            # copy to workspace
            dA[i,j] = sum(view(dA_total_proj, idx1))
            mean_weight_v_no_cb[i,j] = mean(view(projected_velocities_no_cb, idx3))
            mean_weight_v_earth_orb[i,j] = mean(view(v_earth_orb_proj, idx3))

            # calc limb darkening
            for l in eachindex(wavelength)
                ld_sub = map(x -> quad_limb_darkening_eclipse(x, wavelength[l]), mu_grid)
                ld[i,j,l] = sum(view(ld_sub, idx3)) / sum(idx1)
                mean_intensity[i,j,l] = sum(view(ld_sub, idx3)) / sum(idx1)

                # #extinction
                # zenith_angle_matrix = rad2deg.(map(x -> calc_proj_dist(x[1:3], EO_bary[1:3]), OP_bary))
                # extin = map(x -> 10^(-((1/cosd(x))*ext_coef[i,j,l])/2.5), zenith_angle_matrix)
                # ext[i,j,l] = mean(view(extin, idx3)) 

                z_rot[i,j,l] = sum(view(z_rot_sub .* ld_sub .* dA_total_proj, idx3)) ./ sum(view(ld_sub .* dA_total_proj, idx3)) #.* extin

                if isnan(ld[i,j,l])
                    ld[i,j,l] = 0.0
                    ext[i,j,l] = 0.0
                    mean_intensity[i,j,l] = 0.0
                    mean_weight_v_no_cb[i,j] = 0.0
                    mean_weight_v_earth_orb[i,j] = 0.0
                    z_rot[i,j,l] = 0.0
                end
            end
        end
    end
    return 
end

function generate_tloop_eclipse!(tloop::AA{Int}, wsp::SynthWorkspaceEclipse{T}, soldata::SolarData{T}) where T<:AF
    generate_tloop_eclipse!(tloop, wsp.μs, wsp.keys, soldata.len)
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

function get_keys_and_cbs_eclispe!(wsp::SynthWorkspaceEclipse{T}, soldata::SolarData{T}) where T<:AF
    get_keys_and_cbs_eclispe!(wsp.keys, wsp.μs, wsp.cbs, wsp.ax_codes, soldata)
    return nothing
end


function get_keys_and_cbs_eclispe!(keys::AA{Tuple{Symbol, Symbol}}, μs::AA{T1}, cbs::AA{T1},
                           ax_codes::AA{Int}, soldata::SolarData{T2}) where  {T1<:AF, T2<:AF}
    # get the mu and axis codes
    disc_mu = soldata.mu
    disc_ax = soldata.ax

    # type finagling
    disc_mu = convert.(T1, disc_mu)

    # loop over μ positions
    for i in eachindex(μs)
        # move on if we are off the grid
        μs[i] <= zero(T1) && continue

        # get input data for place on disk
        the_key = get_key_for_pos(μs[i], ax_codes[i], disc_mu, disc_ax)
        cbs[i] = convert(T1, soldata.cbs[the_key])
        keys[i] = the_key
    end
    return nothing
end
