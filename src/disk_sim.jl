# line loop function, update prof in place
function line_loop_cpu(prof::AA{T,1}, λΔD::T, depth::T, lambdas::AA{T,1},
                       wsp::SynthWorkspace{T}) where T<:AF
    # first trim the bisectors to the correct depth
    trim_bisector!(depth, wsp.bist, wsp.intt)

    # update the line profile in place
    line_profile_cpu!(λΔD, lambdas, prof, wsp)
    return nothing
end

function time_loop_cpu(tloop::Int, prof::AA{T,1}, z_rot::T, z_cbs::T,
                       z_cbs_avg::T, key::Tuple{Symbol, Symbol},
                       liter::UnitRange{Int}, spec::SpecParams{T},
                       soldata::SolarData, wsp::SynthWorkspace{T}) where T<:AF
    # reset prof
    prof .= one(T)

    # loop over lines
    for l in liter
        # get views needed for line synthesis
        wsp.bist .= copy(view(soldata.bis[key], :, tloop))
        wsp.intt .= copy(view(soldata.int[key], :, tloop))
        wsp.widt .= copy(view(soldata.wid[key], :, tloop))

        # calculate the position of the line center
        extra_z = spec.conv_blueshifts[l] - z_cbs_avg
        λΔD = spec.lines[l] * (1.0 + z_rot) * (1.0 + z_cbs .* spec.variability[l]) * (1.0 + extra_z)

        # get rid of bisector and fix width if variability is turned off
        wsp.bist .*= spec.variability[l]
        if !spec.variability[l]
            wsp.widt .= view(soldata.wid[key], :, 1)
        end

        # get depth to trim to from depth contrast
        dtrim = spec.depths[l] * soldata.dep_contrast[key]

        # synthesize the line
        line_loop_cpu(prof, λΔD, dtrim, spec.lambdas, wsp)
    end
    return nothing
end

function disk_sim(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                  wsp::SynthWorkspace, prof::AA{T,1}, outspec::AA{T,2},
                  tloop::AA{Int,2}; verbose::Bool=true,
                  skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    liter = 1:length(spec.lines); @assert length(liter) >= 1

    # get sum of weights
    sum_wts = sum(wsp.wts)
    z_cbs_avg = sum(wsp.wts .* wsp.cbs) / sum_wts

    # loop over grid position
    for i in CartesianIndices(wsp.μs)
        # move to next iteration if patch element is not visible
        μc = wsp.μs[i]
        μc <= zero(T) && continue

        # get input data for place on disk
        key = wsp.keys[i]
        len = soldata.len[key]

        # get total desired convective blueshift for line
        z_cbs = wsp.cbs[i]

        # get rotational shift
        z_rot = wsp.z_rot[i]

        # loop over time
        for t in 1:disk.Nt
            # check that tloop hasn't exceeded number of epochs
            if tloop[i] > len
                tloop[i] = 1
            end

            # if skip times is true, continue to next iter
            if skip_times[t]
                tloop[i] += 1
                continue
            end

            # update flux profile in place
            time_loop_cpu(tloop[i], prof, z_rot, z_cbs, z_cbs_avg,
                          key, liter, spec, soldata, wsp)

            # apply normalization term and add to outspec
            outspec[:,t] .+= (prof .* wsp.wts[i])

            # iterate tloop
            tloop[i] += 1
        end
    end

    # divide by sum of weights
    outspec ./= sum_wts

    # set instances of outspec where skip is true to 0 and return
    outspec[:, skip_times] .= zero(T)
    return nothing
end
