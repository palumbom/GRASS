function disk_sim_eclipse(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, soldata::SolarData{T},
                            wsp::SynthWorkspaceEclipse{T}, prof::AA{T,1}, flux::AA{T,2},
                            tloop, tloop_init, templates, idx, LD_type, wavelength, 
                            zenith_mean, dA_total_proj, idx1, idx3, mu_grid, z_rot_sub,
                            stored_μs, stored_ax_codes, stored_dA, neid_ext_coeff, ext_toggle; verbose::Bool=true,
                            skip_times::BitVector=falses(disk.Nt)) where T<:AF
    extinction_coeff = DataFrame(CSV.File("data/NEID_two_extinctions.csv"))
    # loop over time
    for t in 1:disk.Nt
            if t < 85
                neid_ext_coeff_new = extinction_coeff[extinction_coeff[!, "Wavelength"] .== wavelength, "Ext1"]
            elseif t >= 85
                neid_ext_coeff_new = extinction_coeff[extinction_coeff[!, "Wavelength"] .== wavelength, "Ext2"]
            end

            #compute intensity for timestamp
            GRASS.eclipse_compute_intensity(disk, wavelength, neid_ext_coeff_new, LD_type, idx1[t], idx3[t],
                                mu_grid[t], z_rot_sub[t], dA_total_proj[t], wsp.ld, wsp.z_rot, zenith_mean[t], 
                                stored_μs, stored_ax_codes, stored_dA, wsp.μs, wsp.ax_codes, wsp.dA, ext_toggle, t, wsp.ext)
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