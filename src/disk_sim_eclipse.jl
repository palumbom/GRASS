function disk_sim_eclipse(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
    wsp::SynthWorkspace{T}, prof::AA{T,1}, flux::AA{T,2},
    tloop::AA{Int,1}; verbose::Bool=true,
    skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # get sum of weights
    sum_wts = sum(wsp.wts)
    z_cbs_avg = sum(wsp.wts .* wsp.cbs) / sum_wts

    # loop over time
    for t in 1:disk.Nt
    # if skip times is true, continue to next iter
    if skip_times[t]
        tloop .+= 1
    continue
    end

    #1: disk.O⃗ function to get this updated for each time step 
    eclipse_compute_quantities(disk, ϕc, θc, μs, ld, dA, xyz, wts, z_rot, ax_codes)

    # loop over wavelength
    for l in 1:length(spec.lines)
    # reset prof
    prof .= zero(T)

    # loop over spatial patches
    for i in eachindex(wsp.μs)
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

    # check that tloop hasn't exceeded number of epochs
    if tloop[i] > len
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
    line_profile_cpu!(λΔD, wsp.wts[i], spec.lambdas, prof, wsp)
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