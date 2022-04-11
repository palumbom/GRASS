using Peaks
using LsqFit
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()

function fit_line_wings(wavs, spec; center=NaN, side="left")
    @assert !isnan(center)

    # find flux at 25 and 80 percent depth
    flux_lo = 0.10 * (1.0 - minimum(spec)) + minimum(spec)
    flux_hi = 0.9 #0.80 * (1.0 - minimum(spec)) + minimum(spec)

    # cut out the middle
    botind = argmin(spec)
    if side == "left"
        idx_lo = findfirst(x -> x .<= flux_lo, spec)
        idx_hi = 1#findfirst(x -> x .<= flux_hi, spec)
        wavs_fit = wavs[idx_hi:idx_lo]
        spec_fit = spec[idx_hi:idx_lo]

        # plt.plot(wavs, spec);
        # plt.axvline(wavs[botind], c="k")
        # plt.axhline(spec[idx_lo], c="k")
        # plt.axhline(spec[idx_hi], c="k")
        # plt.show()
    elseif side == "right"
        idx_lo = findfirst(x -> x .>= flux_lo, spec[botind:end]) + botind
        idx_hi = length(spec)#findfirst(x -> x .>= flux_hi, spec[botind:end]) + botind
        idx_hi = length(spec)
        wavs_fit = wavs[idx_lo:idx_hi]
        spec_fit = spec[idx_lo:idx_hi]

        # plt.clf();
        # plt.plot(wavs, spec);
        # plt.axvline(wavs[botind], c="k")
        # plt.axhline(spec[idx_lo], c="k")
        # plt.axhline(spec[idx_hi], c="k")
    end

    # perform the fit
    p0 = [-0.125, center, 0.03, 0.03, 1.0]
    fit = curve_fit(fit_voigt, wavs_fit, spec_fit, p0)
    return fit_voigt(wavs, fit.param)
end

function clean_line(wavs::AA{T,1}, spec::AA{T,1}; center::T=NaN, plot=false,
                    lwing_Δ_idx::Int=0, rwing_Δ_idx::Int=0) where T<:Float64
    # check that the length of arrays match
    @assert !isnan(center)
    @assert length(wavs) == length(spec)

    # find peaks in spectrum
    m_inds = Peaks.argminima(spec, 12, strict=false)
    m_inds, m_proms = Peaks.peakproms(m_inds, spec, minprom=0.1*std(spec), strict=false)
    m_inds, m_widths, m_left, m_right = Peaks.peakwidths(m_inds, spec, m_proms, strict=false)

    # get better center estimate
    m_inds_idx = argmin(abs.(wavs[m_inds] .- center))
    center_idx = m_inds[m_inds_idx]
    center_wav = wavs[center_idx]

    # find indices for continuum on either side of line
    if iszero(lwing_Δ_idx) && iszero(rwing_Δ_idx)
        c_inds = Peaks.argmaxima(spec, 10, strict=false)
        c_idx_r = findfirst(wavs[c_inds] .> center_wav)
        c_idx_l = c_idx_r - 1

        @show c_inds[c_idx_l]
        @show c_inds[c_idx_r]

        plt.axvline(wavs[c_inds][c_idx_l], c="k")
        plt.axvline(wavs[c_inds][c_idx_r], c="k")
        plt.axvline(center_wav, c="orange")
        plt.plot(wavs, spec./maximum(spec))
        plt.show()
        return nothing
    end

    # calculate wing indices
    lwing_idx = center_idx - lwing_Δ_idx
    rwing_idx = center_idx + rwing_Δ_idx

    if plot
        plt.axvline(wavs[lwing_idx], c="k")
        plt.axvline(wavs[rwing_idx], c="k")
        plt.axvline(center_wav, ls="--", c="k")
        plt.plot(wavs, spec./maximum(spec[lwing_idx:rwing_idx]))
    end

    # take view of spectrum isolated on line to fit
    newwavs = view(wavs, lwing_idx:rwing_idx)
    newspec = view(spec, lwing_idx:rwing_idx)
    newspec ./= maximum(newspec)
    spec ./= maximum(newspec)

    # find fluxes above 90% to replace with fit
    topint = 0.9
    botind = argmin(newspec)
    ind1 = findfirst(x -> x .<= topint, newspec[1:botind])
    ind2 = findfirst(x -> x .>= topint, newspec[botind:end]) + botind

    # abort if indices are weird
    if isnothing(ind1) || isnothing(ind2)
        println("Indices are weird!!")
        wavs .= NaN
        spec .= NaN
        return wavs, spec
    end

    # get model for line wings
    lwing_flux = fit_line_wings(newwavs, newspec, center=center, side="left")
    rwing_flux = fit_line_wings(newwavs, newspec, center=center, side="right")

    if plot
        plt.plot(newwavs[1:ind1], lwing_flux[1:ind1], c="tab:orange")
        plt.plot(newwavs[ind2:end], rwing_flux[ind2:end], c="tab:orange")
        plt.show()
    end

    # replace data wings with model wings
    newspec[1:ind1] .= lwing_flux[1:ind1]
    newspec[ind2:end] .= rwing_flux[ind2:end]

    # divide out slope of spectrum across line to ensure normalization
    slope = (newspec[end] - newspec[1]) / (newwavs[end] - newwavs[1])
    vals = slope .* (newwavs .- newwavs[1]) .+ newspec[1]
    newspec ./= vals

    # now just set rest of spectrum to 1
    spec[1:lwing_idx] .= 1.0
    spec[rwing_idx:end] .= 1.0
    if plot; plt.plot(wavs, spec); plt.show(); end;
    return wavs, spec
end

function clean_line(wavs::AA{T,2}, spec::AA{T,2}; kwargs...) where T<:Float64
    f = (x,y) -> clean_line(x, y; kwargs...)
    out = map(f, eachcol(wavs), eachcol(spec))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end
