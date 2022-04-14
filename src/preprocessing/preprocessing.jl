using Peaks
using LsqFit
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()

function fit_line_wings(wavs, spec; center=NaN, side="left")
    @assert !isnan(center)

    # find flux at percent depths
    botind = argmin(spec)
    flux_lo = 0.25 * (1.0 - minimum(spec)) + minimum(spec)
    flux_hi = 0.80 * (1.0 - minimum(spec)) + minimum(spec)

    plt.axhline(flux_lo, c="k")
    plt.axhline(flux_hi, c="k")

    if side == "left"
        idx_lo = findfirst(x -> x .<= flux_lo, spec)
        idx_hi = findfirst(x -> x .<= flux_hi, spec)
        wavs_fit = wavs[idx_hi:idx_lo]
        spec_fit = spec[idx_hi:idx_lo]

        # plt.plot(wavs, spec);
        # plt.axvline(wavs[botind], c="k")
        # plt.axhline(spec[idx_lo], c="k")
        # plt.axhline(spec[idx_hi], c="k")
        # plt.show()
    elseif side == "right"
        idx_lo = findfirst(x -> x .>= flux_lo, view(spec, botind:length(spec))) + botind
        idx_hi = findfirst(x -> x .>= flux_hi, view(spec, botind:length(spec))) + botind
        wavs_fit = wavs[idx_lo:idx_hi]
        spec_fit = spec[idx_lo:idx_hi]

        # plt.clf();
        # plt.plot(wavs, spec);
        # plt.axvline(wavs[botind], c="k")
        # plt.axhline(spec[idx_lo], c="k")
        # plt.axhline(spec[idx_hi], c="k")
    end

    plt.plot(wavs_fit, spec_fit, c="tab:green")

    # perform the fit
    lb = [-1.0, center-0.05, 0.0, 0.0, 1.0]
    ub = [0.0, center+0.05, 5.0, 5.0, 1.0]
    p0 = [-0.125, center, 0.03, 0.03, 1.0]
    fit = curve_fit(fit_voigt, wavs_fit, spec_fit, p0)# lower=lb, upper=ub)
    # @show fit.param
    return fit_voigt(wavs, fit.param)
end

function clean_line(wavs::AA{T,1}, spec::AA{T,1}; center::T=NaN, plot=false,
                    lwing_Δ_idx::Int=0, rwing_Δ_idx::Int=0) where T<:Float64
    # check that the length of arrays match
    @assert !isnan(center)
    @assert length(wavs) == length(spec)

    # compute a moving average
    n_avg = 5
    wavs_ma = moving_average(wavs, 5)
    spec_ma = moving_average(spec, 5)

    # strip plateaus
    plat_inds = Array{Bool,1}(undef, length(spec_ma))
    plat_inds .= false
    for i in 2:length(spec_ma)
        if isequal(spec_ma[i-1], spec_ma[i])
            plat_inds[i] = true
        end
    end

    # find peaks in spectrum
    m_inds = Peaks.argminima(spec_ma[.!plat_inds], 10, strict=false)

    # get better center estimate
    m_inds_idx = argmin(abs.(wavs[m_inds] .- center))
    center_idx = m_inds[m_inds_idx]
    center_wav = wavs[center_idx]

    # find indices for continuum on either side of line
    if iszero(lwing_Δ_idx) && iszero(rwing_Δ_idx)
        c_inds = Peaks.argmaxima(spec, 10, strict=false)
        c_idx_r = findfirst(wavs[c_inds] .> center_wav)
        c_idx_l = c_idx_r - 1

        @show center_idx - c_inds[c_idx_l]
        @show c_inds[c_idx_r] - center_idx

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

    # do normalization
    spec ./= maximum(newspec)
    newspec ./= maximum(newspec)

    # # find fluxes above 90% to replace
    # topint = 0.9
    # botind = argmin(newspec);
    # # @show botind
    # # if botind == 134
    # #     for i in m_inds
    # #         plt.axvline(wavs[i])
    # #     end
    # #     plt.plot(wavs_ma, spec_ma)
    # #     plt.plot(newwavs, newspec)
    # #     plt.show()
    # # end
    # ind1 = findfirst(x -> x .<= topint, newspec[1:botind])
    # ind2 = findfirst(x -> x .>= topint, newspec[botind:end]) + botind

    # # abort if indices are weird
    # if isnothing(ind1) || isnothing(ind2)
    #     println("Indices are weird!!")
    #     wavs .= NaN
    #     spec .= NaN
    #     return wavs, spec
    # end

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
