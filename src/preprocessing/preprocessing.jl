using LsqFit
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()

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

function isolate_line(wavs::AA{T,1}, spec::AA{T,1}; center::T=NaN) where T<:Float64
    # check that the length of arrays match
    @assert !isnan(center)
    @assert length(wavs) == length(spec)

    # get better center estimate
    if (center - 1.0) < minimum(wavs)
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

    # divide out slope of spectrum across line to ensure normalization
    slope = (newspec[end] - newspec[1]) / (newwavs[end] - newwavs[1])
    vals = slope .* (newwavs .- newwavs[1]) .+ newspec[1]
    newspec ./= vals
    return newwavs, newspec
end

function isolate_line(wavs::AA{T,2}, spec::AA{T,2}; kwargs...) where T<:Float64
    f = (x,y) -> isolate_line(x, y; kwargs...)
    out = map(f, eachcol(wavs), eachcol(spec))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end
