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
        wavs_fit = wavs[botind+15:end-15]
        spec_fit = spec[botind+15:end-15]
    end

    # perform the fit
    p0 = [-0.125, center, 0.03, 0.03, 1.0]
    fit = curve_fit(fit_voigt, wavs_fit, spec_fit, p0)
    return fit_voigt(wavs, fit.param)
end

function clean_line(wavs::AA{T,1}, spec::AA{T,1}; center::T=NaN, plot=false) where T<:Float64
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
    buffer = 1.25
    lwingλ = center - buffer
    rwingλ = center + buffer

    # indices to narrow region of interest
    lwavind = searchsortednearest(wavs, lwingλ)
    rwavind = searchsortednearest(wavs, rwingλ)

    # check the inds
    if isnothing(lwavind)
        lwavind = 1
    end

    if isnothing(rwavind)
        rwavind = length(wavs)
    end

    if plot; plt.plot(wavs, spec./maximum(spec[lwavind:rwavind])); end

    # take view of spectrum isolated on line
    newwavs = view(wavs, lwavind:rwavind)
    newspec = view(spec, lwavind:rwavind)
    newspec ./= maximum(newspec)

    # find the wings
    topint = 0.9
    botind = argmin(newspec)
    ind1 = findfirst(x -> x .<= topint, newspec[1:botind])
    ind2 = findfirst(x -> x .>= topint, newspec[botind:end]) + botind

    # abort if indices are weird
    if isnothing(ind1) || isnothing(ind2)
        wavs .= NaN
        spec .= NaN
        return wavs, spec
    end

    # get model for line wings
    lwing_flux = fit_line_wings(newwavs, newspec, center=center, side="left")
    rwing_flux = fit_line_wings(newwavs, newspec, center=center, side="right")

    # replace data wings with model wings
    newspec[1:ind1] .= lwing_flux[1:ind1]
    newspec[ind2:end] .= rwing_flux[ind2:end]

    # divide out slope of spectrum across line to ensure normalization
    slope = (newspec[end] - newspec[1]) / (newwavs[end] - newwavs[1])
    vals = slope .* (newwavs .- newwavs[1]) .+ newspec[1]
    newspec ./= vals

    # normalize it one last time
    newspec ./= maximum(newspec)

    # now just set rest of spectrum to 1
    spec[1:lwavind] .= 1.0
    spec[rwavind:end] .= 1.0
    if plot; plt.plot(wavs, spec); plt.show(); end;
    return wavs, spec
end

function clean_line(wavs::AA{T,2}, spec::AA{T,2}; kwargs...) where T<:Float64
    f = (x,y) -> clean_line(x, y; kwargs...)
    out = map(f, eachcol(wavs), eachcol(spec))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end
