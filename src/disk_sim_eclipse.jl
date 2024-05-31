function disk_sim_eclipse(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, soldata::SolarData{T},
    wsp::SynthWorkspaceEclipse{T}, mem::GeoWorkspaceEclipse{T}, prof::AA{T,1}, flux::AA{T,2},
    tloop, tloop_init, templates, idx, obs_long, obs_lat, alt, time_stamps, band, wavelength, CB::Bool; verbose::Bool=true,
    skip_times::BitVector=falses(disk.Nt)) where T<:AF

    # loop over time
    for t in 1:disk.Nt

            #compute geometry for timestamp
            GRASS.eclipse_compute_quantities!(disk, time_stamps[t], band, obs_long, obs_lat, alt, wavelength, wsp.ϕc, wsp.θc, 
                                                wsp.μs, wsp.ld, wsp.ext, wsp.dA, wsp.xyz, wsp.z_rot, wsp.ax_codes,
                                                mem.dA_total_proj_mean, mem.mean_intensity, mem.mean_weight_v_no_cb,
                                                mem.mean_weight_v_earth_orb, mem.pole_vector_grid,
                                                mem.SP_sun_pos, mem.SP_sun_vel, mem.SP_bary, mem.SP_bary_pos,
                                                mem.SP_bary_vel, mem.OP_bary, mem.mu_grid, mem.projected_velocities_no_cb, 
                                                mem.distance, mem.v_scalar_grid, mem.v_earth_orb_proj)
            # get conv. blueshift and keys from input data
            GRASS.get_keys_and_cbs_eclispe!(wsp, soldata)
            
            # generate or copy tloop
            if (idx > 1) && GRASS.in_same_group(templates[idx - 1], templates[idx])
                tloop .= tloop_init
            else
                GRASS.generate_tloop_eclipse!(tloop_init, wsp, soldata)
                tloop .= tloop_init
            end

            # loop over wavelength
            for l in 1:length(spec.lines)
                # reset prof
                prof .= zero(T)

                # get sum of weights
                sum_wts = sum(wsp.ld[:, :, l] .* wsp.dA) #.*wsp.ext
                z_cbs_avg = sum(wsp.ld[:, :, l] .* wsp.dA .* wsp.cbs) / sum_wts

                contrast = (wsp.ld[:, :, l] / NaNMath.maximum(wsp.ld[:, :, l])).^0.1

                # loop over spatial patches
                for i in eachindex(disk.ϕc)
                    for j in 1:disk.Nθ[i]

                    # move to next iteration if patch element is not visible
                    (wsp.ld[i,j,l] .* wsp.dA[i,j]) <= zero(T) && continue

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

                    if CB == true
                        # get amount of convective blueshift needed
                        extra_z = spec.conv_blueshifts[l] - z_cbs_avg
                    else
                        extra_z = 0.0
                    end

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
                    line_profile_cpu!(λΔD, wsp.dA[i,j], wsp.ld[i,j,l], wsp.ext[i,j,l], contrast[i,j], spec.lambdas, prof, wsp)
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