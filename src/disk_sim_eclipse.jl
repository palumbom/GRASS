function disk_sim_eclipse(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, soldata::SolarData{T},
                            wsp::SynthWorkspaceEclipse{T}, prof::AA{T,1}, flux::AA{T,2},
                            tloop, tloop_init, templates::AA{String,1}, idx, LD_type::String, wavelength::Vector{Float64}, 
                            time_stamps::Vector{Float64}, obs_long::T, obs_lat::T, alt::T,
                            ext_coeff, ext_toggle::Bool; skip_times::BitVector=falses(disk.Nt)) where T<:AF

    # loop over time
    for t in 1:disk.Nt
            GRASS.eclipse_compute_quantities(time_stamps[t], t, obs_long, obs_lat, alt, wavelength, LD_type, ext_toggle, ext_coeff, disk, wsp)

            # get conv. blueshift and keys from input data
            GRASS.get_keys_and_cbs_eclispe!(wsp, soldata, t)
            
            # generate or copy tloop
            if (idx > 1) && GRASS.in_same_group(templates[idx - 1], templates[idx])
                tloop .= tloop_init
            else
                GRASS.generate_tloop_eclipse!(tloop_init, wsp, soldata, t)
                tloop .= tloop_init
            end

            # loop over wavelength
            for l in 1:length(spec.lines)
                # reset prof
                prof .= zero(T)

                # get sum of weights
                if ext_toggle == false
                    sum_wts = sum(wsp.ld[:, :, l] .* wsp.dA[:, :, t])
                    z_cbs_avg = sum(wsp.ld[:, :, l] .* wsp.dA[:, :, t] .* wsp.cbs) / sum_wts
                end
                if ext_toggle == true
                    sum_wts = sum(wsp.ld[:, :, l] .* wsp.dA[:, :, t] .* wsp.ext[:, :, l])
                    z_cbs_avg = sum(wsp.ld[:, :, l] .* wsp.dA[:, :, t] .* wsp.ext[:, :, l] .* wsp.cbs) / sum_wts
                end

                # loop over spatial patches
                for i in eachindex(disk.ϕc)
                    for j in 1:disk.Nθ[i]

                    # move to next iteration if patch element is not visible
                    if ext_toggle == false
                        (wsp.ld[i,j,l] .* wsp.dA[i,j,t]) <= zero(T) && continue
                    end
                    if ext_toggle == true
                        (wsp.ld[i,j,l] .* wsp.dA[i,j,t] .* wsp.ext[i,j,l]) <= zero(T) && continue
                    end

                    # get input data for place on disk
                    key = wsp.keys[i,j]
                    len = soldata.len[key]

                    # get total desired convective blueshift for line
                    z_cbs = wsp.cbs[i,j]

                    # get rotational shift
                    z_rot = wsp.z_rot[i,j,l]

                    # check that tloop hasn't exceeded number of epochs
                    if tloop[i,j] > len
                        tloop[i,j] -= len
                    end

                    # get views needed for line synthesis
                    wsp.bist .= copy(view(soldata.bis[key], :, tloop[i,j]))
                    wsp.intt .= copy(view(soldata.int[key], :, tloop[i,j])) 
                    wsp.widt .= copy(view(soldata.wid[key], :, tloop[i,j]))

                    # get amount of convective blueshift needed
                    extra_z = spec.conv_blueshifts[l] - z_cbs_avg

                    # get shifted line center
                    λΔD = spec.lines[l]
                    λΔD *= (1.0 + z_rot)
                    λΔD *= (1.0 + z_cbs .* spec.variability[l])
                    λΔD *= (1.0 + extra_z .* spec.variability[l])

                    # get rid of bisector and fix width if variability is turned off
                    wsp.bist .*= spec.variability[l]
                    if !spec.variability[l]
                        wsp.widt .= view(soldata.wid[key], :, 1)
                    end

                    # get depth to trim to from depth contrast
                    dtrim = spec.depths[l] * soldata.dep_contrast[key]

                    # first trim the bisectors to the correct depth
                    trim_bisector!(dtrim, wsp.bist, wsp.intt)

                    # update the line profile in place
                    line_profile_cpu!(λΔD, wsp.dA[i,j,t], wsp.ld[i,j,l], wsp.ext[i,j,l], spec.lambdas, prof, wsp, ext_toggle)
                    end
                end

                # apply normalization term and add to flux
                flux[:,t] .*= prof./ sum_wts 
            end

            # iterate tloop
            tloop .+= 1
    end

    # set instances of outspec where skip is true to 0 and return
    flux[:, skip_times] .= zero(T)
    return nothing
end

function disk_sim_rossiter(spec::SpecParams{T}, disk::DiskParams{T}, planet::Planet{T},
                           soldata::SolarData{T}, wsp::SynthWorkspace{T}, ros_allocs::RossiterAllocs{T},
                           prof::AA{T,1}, flux::AA{T,2}, vels::AA{T,1}, tloop::AA{Int,1}; verbose::Bool=true,
                           skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # get sum of weights
    sum_wts_og = sum(wsp.wts)
    z_cbs_avg = sum(wsp.wts .* wsp.cbs) / sum_wts_og

    # alias allocations
    xyz_planet = ros_allocs.xyz_planet
    xyz_dot_planet = ros_allocs.xyz_dot_planet
    xyz_star = ros_allocs.xyz_star
    xyz_dot_star = ros_allocs.xyz_dot_star

    # loop over time
    for t in 1:disk.Nt
        # if skip times is true, continue to next iter
        if skip_times[t]
            tloop .+= 1
            continue
        end

        # recopy weights, etc. from unobstructed disk
        ros_allocs.μs .= wsp.μs
        ros_allocs.ld .= wsp.ld
        ros_allocs.dA .= wsp.dA
        ros_allocs.wts .= wsp.wts
        ros_allocs.z_rot .= wsp.z_rot

        # calculate projected distance between planet and star centers
        dist2 = calc_proj_dist2(xyz_planet[:,t], xyz_star[:,t])

        # check if the planet is transiting
        if (dist2 <= (disk.ρs + planet.radius)^2.0) & (xyz_planet[3] > 0.0)
            # re-calculate patch weights, etc. for occulted patches
            calc_rossiter_quantities!(xyz_planet[:,t], planet, disk, wsp, ros_allocs)
        end

        # get the sum of the weights w/ occulted pixels
        sum_wts = sum(ros_allocs.wts)

        # get weighted sum of velocities
        vels[t] = sum(ros_allocs.z_rot .* c_ms .* ros_allocs.wts) ./ sum_wts

        # loop over wavelength
        for l in 1:length(spec.lines)
            # reset prof
            prof .= zero(T)

            # loop over spatial patches
            for i in eachindex(wsp.μs)
                # move to next patch if entire patch is behind hemisphere
                μc = ros_allocs.μs[i]
                μc <= zero(T) && continue

                # move to next iteration if entire patch is occulted
                wts = ros_allocs.wts[i]
                iszero(wts) && continue

                # get input data for place on disk
                key = wsp.keys[i]
                len = soldata.len[key]

                # get total desired convective blueshift for line
                z_cbs = wsp.cbs[i]

                # get rotational shift
                z_rot = ros_allocs.z_rot[i]

                # check that tloop hasn't exceeded number of epochs
                while tloop[i] > len
                    tloop[i] -= len
                end

                # get views needed for line synthesis
                wsp.bist .= copy(view(soldata.bis[key], :, tloop[i]))
                wsp.intt .= copy(view(soldata.int[key], :, tloop[i]))
                wsp.widt .= copy(view(soldata.wid[key], :, tloop[i]))

                # get amount of convective blueshift needed
                extra_z = spec.conv_blueshifts[l] - z_cbs_avg

                # get shifted line center
                λΔD = spec.lines[l]
                λΔD *= (1.0 + z_rot)
                λΔD *= (1.0 + z_cbs .* spec.variability[l])
                λΔD *= (1.0 + extra_z)

                # get rid of bisector and fix width if variability is turned off
                wsp.bist .*= spec.variability[l]
                if !spec.variability[l]
                    wsp.widt .= view(soldata.wid[key], :, 1)
                end

                # get depth to trim to from depth contrast
                dtrim = spec.depths[l] * soldata.dep_contrast[key]

                # first trim the bisectors to the correct depth
                trim_bisector!(dtrim, wsp.bist, wsp.intt)

                # update the line profile in place
                line_profile_cpu!(λΔD, wts, spec.lambdas, prof, wsp)
            end

            # apply normalization term and add to flux
            flux[:,t] .*= prof ./ sum_wts
        end

        # iterate tloop
        tloop .+= 1
    end

    # set instances of outspec where skip is true to 0 and return
    flux[:, skip_times] .= zero(T)
    return nothing
end
