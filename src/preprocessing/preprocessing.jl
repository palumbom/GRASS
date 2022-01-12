using LsqFit

function calc_width_at_depth(wavs::AA{T,1}, spec::AA{T,1}; center::T=NaN,
                             len::Integer=100) where T<:Real
    # check that the length of arrays match
    @assert !isnan(center)
    @assert length(wavs) == length(spec)

    # get better center estimate
    if (center - 1.0) > minimum(wavs)
        lwingλ_idx = firstindex(wavs)
    else
        lwingλ_idx = findfirst(x -> x >= center - 1.0, wavs)
    end

    if (center + 1.0) > maximum(wavs)
        rwingλ_idx = lastindex(wavs)
    else
        rwingλ_idx = findfirst(x -> x >= center + 1.0, wavs)
    end
    center_idx = argmin(view(spec, lwingλ_idx:rwingλ_idx)) + lwingλ_idx
    center = wavs[center_idx]

    # get better wing estimate
    buffer = 0.4
    lwingλ = center - buffer
    rwingλ = center + buffer

    # hard-coded indices to narrow region of interest
    lwavind = searchsortednearest(wavs, lwingλ)
    rwavind = searchsortednearest(wavs, rwingλ)
    newwavs = wavs[lwavind:rwavind]
    newspec = spec[lwavind:rwavind]
    newspec ./= maximum(newspec)

    # find the wings
    topint = 0.93
    botind = argmin(newspec)
    ind1 = findfirst(x -> x .<= topint, newspec[1:botind])
    ind2 = findfirst(x -> x .>= topint, newspec[botind:end]) + botind

    # get model for line wings
    lwing_wavs, lwing_flux = fit_line_wings(newwavs, newspec, center=center, side="left")
    rwing_wavs, rwing_flux = fit_line_wings(newwavs, newspec, center=center, side="right")

    # replace data wings with model wings
    newwavs[1:ind1] .= lwing_wavs[1:ind1]
    newspec[1:ind1] .= lwing_flux[1:ind1]
    newwavs[ind2:end] .= rwing_wavs[ind2:end]
    newspec[ind2:end] .= rwing_flux[ind2:end]

    # make sure it's normalized one more time
    newspec ./= maximum(newspec)

    # measure the width via interpolation method and return
    dep, wid = measure_width_interpolate(newwavs, newspec)
    return dep, wid
end

function calc_width_at_depth(wavs::AA{T,2}, spec::AA{T,2}; center::T=NaN,
                             len::Integer=100) where T<:Real
    f = (x,y) -> calc_width_at_depth(x, y, len=len, center=center)
    out = map(f, eachcol(wavs), eachcol(spec))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end

function fit_line_wings(wavs, spec; center=NaN, side="left")
    @assert !isnan(center)

    # cut out the middle
    botind = argmin(spec)
    if side == "left"
        wavs_fit = wavs[15:botind-15]
        spec_fit = spec[15:botind-15]
    elseif side == "right"
        wavs_fit = wavs[botind+10:end-15]
        spec_fit = spec[botind+10:end-15]
    end

    # perform the fit
    p0 = [-0.125, center, 0.03, 0.03, 1.0]
    fit = curve_fit(fit_voigt, wavs_fit, spec_fit, p0)
    return wavs, fit_voigt(wavs, fit.param)
end

function measure_width_interpolate(xs::AA{T,1}, ys::AA{T,1}; top::T=1.0,
                                   len::Integer=100) where T<:AF
    # check lengths and normalization
    @assert length(xs) == length(ys)

    # normalize the spec, find bottom of line
    # ys ./= maximum(ys)
    botind = argmin(ys)
    depths = range(ys[botind], top, length=len)

    # find left and right halves
    lind = findfirst(ys .< top)
    rind = findlast(ys .< top)
    lspec = reverse(ys[lind:botind])
    rspec = ys[botind:rind]
    lwav = reverse(xs[lind:botind])
    rwav = xs[botind:rind]

    # interpolate onto even intensity grid and compute width
    lspline = linear_interp(lspec, lwav)
    rspline = linear_interp(rspec, rwav)
    wids = rspline.(depths) .- lspline.(depths)
    return depths, wids
end

function calc_bisector_at_depth(wavs::AA{T,1}, spec::AA{T,1}; center::T=NaN,
                                len::Integer=100) where T<:Real
    # check that the length of arrays match
    @assert !isnan(center)
    @assert length(wavs) == length(spec)

    # get better center estimate
    if (center - 1.0) > minimum(wavs)
        lwingλ_idx = firstindex(wavs)
    else
        lwingλ_idx = findfirst(x -> x >= center - 1.0, wavs)
    end

    if (center + 1.0) > maximum(wavs)
        rwingλ_idx = lastindex(wavs)
    else
        rwingλ_idx = findfirst(x -> x >= center + 1.0, wavs)
    end
    center_idx = argmin(view(spec, lwingλ_idx:rwingλ_idx)) + lwingλ_idx
    center = wavs[center_idx]

    # get better wing estimate
    buffer = 0.4
    lwingλ = center - buffer
    rwingλ = center + buffer

    # hard-coded indices to narrow region of interest
    lwavind = searchsortednearest(wavs, lwingλ)
    rwavind = searchsortednearest(wavs, rwingλ)
    newwavs = wavs[lwavind:rwavind]
    newspec = spec[lwavind:rwavind]
    newspec ./= maximum(newspec)

    # find the wings
    topint = 0.93
    botind = argmin(newspec)
    ind1 = findfirst(x -> x .<= topint, newspec[1:botind])
    ind2 = findfirst(x -> x .>= topint, newspec[botind:end]) + botind

    # get model for line wings
    lwing_wavs, lwing_flux = fit_line_wings(newwavs, newspec, center=center, side="left")
    rwing_wavs, rwing_flux = fit_line_wings(newwavs, newspec, center=center, side="right")

    # replace data wings with model wings
    newwavs[1:ind1] .= lwing_wavs[1:ind1]
    newspec[1:ind1] .= lwing_flux[1:ind1]
    newwavs[ind2:end] .= rwing_wavs[ind2:end]
    newspec[ind2:end] .= rwing_flux[ind2:end]

    # make sure it's normalized one more time
    newspec ./= maximum(newspec)

    # measure the width via interpolation method and return
    wav, bis = measure_bisector_interpolate(newwavs, newspec)
    return wav, bis
end

function calc_bisector_at_depth(wavs::AA{T,2}, spec::AA{T,2}; center::T=NaN,
                                len::Integer=100) where T<:Real
    f = (x,y) -> calc_bisector_at_depth(x, y, len=len, center=center)
    out = map(f, eachcol(wavs), eachcol(spec))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end
