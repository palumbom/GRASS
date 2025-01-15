function disk_sim_eclipse(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, soldata::SolarData{T},
                            wsp::SynthWorkspaceEclipse{T}, prof::AA{T,1}, flux::AA{T,2},
                            tloop, tloop_init, templates::AA{String,1}, idx, LD_type::String, wavelength::Vector{Float64}, 
                            zenith_mean::Vector{Matrix{Matrix{Float64}}}, dA_total_proj::Vector{Matrix{Matrix{Float64}}}, 
                            idx1::Vector{Matrix{Matrix{Int}}}, idx3::Vector{Matrix{Matrix{Int}}}, mu_grid::Vector{Matrix{Matrix{Float64}}}, 
                            z_rot_sub::Vector{Matrix{Matrix{Float64}}}, stored_μs::AA{T,3}, stored_ax_codes::AA{Int,3}, 
                            stored_dA::AA{T,3}, ext_coeff, ext_toggle::Bool;
                            skip_times::BitVector=falses(disk.Nt)) where T<:AF

    # loop over time
    for t in 1:disk.Nt
            if ext_coeff == "three"
                if LD_type == "SSD"
                    ext_coeff_file = ext_file_KSSD
                elseif LD_type == "300"
                    ext_coeff_file = ext_file_K300
                end
                
                if t < 25
                    coeff1 = ext_coeff_file[ext_coeff_file[!, "Wavelength"] .== wavelength, "Ext1"][1]
                    # compute intensity for timestamp
                    GRASS.eclipse_compute_intensity(disk, wavelength, coeff1, LD_type, idx1[t], idx3[t],
                                mu_grid[t], z_rot_sub[t], dA_total_proj[t], wsp.ld, wsp.z_rot, zenith_mean[t], 
                                stored_μs, stored_ax_codes, stored_dA, wsp.μs, wsp.ax_codes, wsp.dA, ext_toggle, t, wsp.ext)
                elseif t >= 25 && t < 46 
                    coeff2 = ext_coeff_file[ext_coeff_file[!, "Wavelength"] .== wavelength, "Ext2"][1]
                    # compute intensity for timestamp
                    GRASS.eclipse_compute_intensity(disk, wavelength, coeff2, LD_type, idx1[t], idx3[t],
                                mu_grid[t], z_rot_sub[t], dA_total_proj[t], wsp.ld, wsp.z_rot, zenith_mean[t], 
                                stored_μs, stored_ax_codes, stored_dA, wsp.μs, wsp.ax_codes, wsp.dA, ext_toggle, t, wsp.ext)
                elseif t >= 46
                    coeff3 = ext_coeff_file[ext_coeff_file[!, "Wavelength"] .== wavelength, "Ext3"][1]
                    # compute intensity for timestamp
                    GRASS.eclipse_compute_intensity(disk, wavelength, coeff3, LD_type, idx1[t], idx3[t],
                                mu_grid[t], z_rot_sub[t], dA_total_proj[t], wsp.ld, wsp.z_rot, zenith_mean[t], 
                                stored_μs, stored_ax_codes, stored_dA, wsp.μs, wsp.ax_codes, wsp.dA, ext_toggle, t, wsp.ext)
                end
            else
                # compute intensity for timestamp
                GRASS.eclipse_compute_intensity(disk, wavelength, ext_coeff, LD_type, idx1[t], idx3[t],
                        mu_grid[t], z_rot_sub[t], dA_total_proj[t], wsp.ld, wsp.z_rot, zenith_mean[t], 
                        stored_μs, stored_ax_codes, stored_dA, wsp.μs, wsp.ax_codes, wsp.dA, ext_toggle, t, wsp.ext)
            end

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
                    λΔD *= (1.0 + z_cbs)
                    λΔD *= (1.0 + extra_z)

                    # fix bisector and width if variability is turned off 
                    if !spec.variability[l]
                        wsp.bist .= mean(soldata.bis[key], dims=2) 
                        wsp.widt .= mean(soldata.wid[key], dims=2)
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