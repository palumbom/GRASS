# Time mean over epochs of each tile's trimmed bisector, intensity grid and width for line l,
# with the same variability gating disk_sim_eclipse applies per cell. The trim puts every tile
# on a common fractional-depth index, which is what lets tiles be blended index by index.
# Summed in epoch order and divided by the count, matching time_average_bis! on the GPU.
function trimmed_time_means(spec::SpecParams{T}, l::Int, soldata::SolarData{T}) where T<:AF
    bis_mean = Dict{Tuple{Symbol,Symbol}, Vector{T}}()
    int_mean = Dict{Tuple{Symbol,Symbol}, Vector{T}}()
    wid_mean = Dict{Tuple{Symbol,Symbol}, Vector{T}}()
    for key in keys(soldata.len)
        len = soldata.len[key]
        ndep = size(soldata.bis[key], 1)
        dtrim = spec.depths[l] * soldata.dep_contrast[key]
        bist = zeros(T, ndep)
        intt = zeros(T, ndep)
        bis_sum = zeros(T, ndep)
        int_sum = zeros(T, ndep)
        wid_sum = zeros(T, ndep)
        for e in 1:len
            bist .= view(soldata.bis[key], :, e) .* spec.variability[l]
            intt .= view(soldata.int[key], :, e)
            GRASS.trim_bisector!(dtrim, bist, intt)
            bis_sum .+= bist
            int_sum .+= intt
            wid_sum .+= view(soldata.wid[key], :, spec.variability[l] ? e : 1)
        end
        bis_mean[key] = bis_sum ./ len
        int_mean[key] = int_sum ./ len
        wid_mean[key] = wid_sum ./ len
    end
    return bis_mean, int_mean, wid_mean
end

# Axis-pooled copy of the time means: every tile at a mu level gets the unweighted mean over
# that level's tiles. Summed in the tile order of the dictionary keys sorted by axis code,
# which is the GPU's tile order, and divided by the count, matching pool_axis_means_gpu!.
function pool_axis_means(means::NTuple{3, Dict{Tuple{Symbol,Symbol}, Vector{T}}}) where T<:AF
    bm, im, wm = means
    ks = sort(collect(keys(bm)), by=k -> (GRASS.parse_mu_string(k[2]), GRASS.parse_ax_string(k[1])))
    bp = Dict{Tuple{Symbol,Symbol}, Vector{T}}()
    ip = Dict{Tuple{Symbol,Symbol}, Vector{T}}()
    wp = Dict{Tuple{Symbol,Symbol}, Vector{T}}()
    for k in ks
        level = [k2 for k2 in ks if k2[2] == k[2]]
        bs = zeros(T, length(bm[k])); is = zeros(T, length(bm[k])); ws = zeros(T, length(bm[k]))
        for k2 in level
            bs .+= bm[k2]; is .+= im[k2]; ws .+= wm[k2]
        end
        bp[k] = bs ./ length(level)
        ip[k] = is ./ length(level)
        wp[k] = ws ./ length(level)
    end
    return bp, ip, wp
end

function disk_sim_eclipse(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, soldata::SolarData{T},
                            wsp::SynthWorkspaceEclipse{T}, prof::AA{T,1}, flux::AA{T,2},
                            tloop, tloop_init, templates::AA{String,1}, idx, LD_type::String, wavelength::Vector{Float64},
                            time_stamps::Vector{String}, obs_long::T, obs_lat::T, alt::T,
                            ext_coeff, ext_toggle::Bool; skip_times::BitVector=falses(disk.Nt),
                            interp_mu::Bool=false, pool_axes::Bool=false) where T<:AF

    # interp_mu / pool_axes: per line, the time-mean trimmed shape of every tile (`means`) and
    # the source the cell is corrected toward (`src_means`: axis-pooled, or the same dicts)
    adjust_means = interp_mu | pool_axes
    means = Vector{NTuple{3, Dict{Tuple{Symbol,Symbol}, Vector{T}}}}()
    src_means = Vector{NTuple{3, Dict{Tuple{Symbol,Symbol}, Vector{T}}}}()
    if adjust_means
        for l in eachindex(spec.lines)
            m = trimmed_time_means(spec, l, soldata)
            push!(means, m)
            push!(src_means, pool_axes ? pool_axis_means(m) : m)
        end
    end

    # loop over time
    for t in 1:disk.Nt
            GRASS.Eclipse.eclipse_compute_quantities!(time_stamps[t], t, obs_long, obs_lat, alt, wavelength, LD_type, ext_toggle, ext_coeff, disk, wsp)

            # get conv. blueshift and keys from input data
            GRASS.Eclipse.get_keys_and_cbs_eclispe!(wsp, soldata, t)

            # interp_mu: bracketing tiles, weights, and the interpolated blueshift
            if interp_mu
                GRASS.Eclipse.get_interp_keys_eclipse!(wsp, soldata, t)
            end

            # generate or copy tloop
            if (idx > 1) && GRASS.in_same_group(templates[idx - 1], templates[idx])
                tloop .= tloop_init
            else
                GRASS.Eclipse.generate_tloop_eclipse!(tloop_init, wsp, soldata, t)
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

                    # wrap tloop into the data length of this tile's current key
                    tloop[i,j] = mod1(tloop[i,j], len)

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
                    GRASS.trim_bisector!(dtrim, wsp.bist, wsp.intt)

                    # interp_mu / pool_axes: replace this tile's time-mean shape by the target
                    # (limb-angle interpolant of the source means, or the source at this
                    # tile), keeping the epoch's departure from the tile's own mean. Same
                    # expression and evaluation order as fill_workspaces_2D_eclipse!.
                    if adjust_means
                        bm, im, wm = means[l]
                        bs, is, ws = src_means[l]
                        key_lo = interp_mu ? wsp.keys_lo[i,j] : key
                        key_hi = interp_mu ? wsp.keys_hi[i,j] : key
                        w = interp_mu ? wsp.wts[i,j] : zero(T)
                        bm_lo = bs[key_lo]; bm_hi = bs[key_hi]; bm_k = bm[key]
                        im_lo = is[key_lo]; im_hi = is[key_hi]; im_k = im[key]
                        wm_lo = ws[key_lo]; wm_hi = ws[key_hi]; wm_k = wm[key]
                        for n in eachindex(wsp.bist)
                            wsp.bist[n] += ((one(T) - w) * bm_lo[n] + w * bm_hi[n]) - bm_k[n]
                            wsp.intt[n] += ((one(T) - w) * im_lo[n] + w * im_hi[n]) - im_k[n]
                            wsp.widt[n] += ((one(T) - w) * wm_lo[n] + w * wm_hi[n]) - wm_k[n]
                        end
                    end

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
