# follows from implementation at https://github.com/ACCarnall/SpectRes/blob/master/spectres/spectral_resampling.py
function get_bin_edges(arr::AA{T,1}) where T<:AF
    edges = zeros(length(arr)+1)
    edges[1] = arr[1] - 0.5 * (arr[2] - arr[1])
    edges[end] = arr[end] + 0.5 * (arr[end] - arr[end-1])
    edges[2:end-1] = 0.5 .* (arr[2:end] .+ arr[1:end-1])
    widths = diff(edges)
    return edges, widths
end

# follows from implementation at https://github.com/ACCarnall/SpectRes/blob/master/spectres/spectral_resampling.py
# TODO algorithm might cause shift based on input wavelength grid?
# might cause issue if wavelengths shift across
function rebin_spectrum(xs_old::AA{T,1}, ys_old::AA{T,1}, xs_new::AA{T,1}) where T<:AF
    # get edges of bins
    old_edges, old_widths = get_bin_edges(xs_old)
    new_edges, new_widths = get_bin_edges(xs_new)

    # allocate memory for output array
    ys_new = zeros(length(xs_new))

    # loop over new bins
    start = 0
    stop = 0
    for i in eachindex(xs_new)
        # boundary conditions
        if new_edges[i] < first(old_edges)
            ys_new[i] = first(ys_old)
            continue
        elseif new_edges[i+1] > last(old_edges)
            ys_new[i] = last(ys_old)
            continue
        end

        while old_edges[start+1] <= new_edges[i]
            start += 1
        end

        while old_edges[stop+1] < new_edges[i+1]
            stop += 1
        end

        if start == stop
            ys_new[i] = ys_old[start]
        else
            start_factor = ((old_edges[start+1] - new_edges[i]) / (old_edges[start+1] - old_edges[start]))
            stop_factor = ((new_edges[i+1] - old_edges[stop]) / (old_edges[stop+1] - old_edges[stop]))

            old_widths[start] *= start_factor
            old_widths[stop] *= stop_factor

            f_widths = old_widths[start:stop] .* ys_old[start:stop]
            ys_new[i] = sum(f_widths)
            ys_new[i] /= sum(old_widths[start:stop])

            old_widths[start] /= start_factor
            old_widths[stop] /= stop_factor
        end
    end
    return ys_new
end

# follows from implementation at https://github.com/ACCarnall/SpectRes/blob/master/spectres/spectral_resampling.py
function rebin_spectrum(xs_old::AA{T,1}, ys_old::AA{T,1}, σs_old::AA{T,1}, xs_new::AA{T,1}, ) where T<:AF
    @assert length(σs_old) == length(ys_old)

    # get edges of bins
    old_edges, old_widths = get_bin_edges(xs_old)
    new_edges, new_widths = get_bin_edges(xs_new)

    # allocate memory for output arrays
    ys_new = zeros(length(xs_new))
    σs_new = zeros(length(xs_new))

    # loop over new bins
    start = 0
    stop = 0
    for i in eachindex(xs_new)
        # boundary conditions
        if new_edges[i] < first(old_edges)
            ys_new[i] = first(ys_old)
            σs_new[i] = first(σs_old)
            continue
        elseif new_edges[i+1] > last(old_edges)
            ys_new[i] = last(ys_old)
            σs_new[i] = last(σs_old)
            continue
        end

        while old_edges[start+1] <= new_edges[i]
            start += 1
        end

        while old_edges[stop+1] < new_edges[i+1]
            stop += 1
        end

        if start == stop
            ys_new[i] = ys_old[start]
            σs_new[i] = σs_old[start]
        else
            start_factor = ((old_edges[start+1] - new_edges[i]) / (old_edges[start+1] - old_edges[start]))
            stop_factor = ((new_edges[i+1] - old_edges[stop]) / (old_edges[stop+1] - old_edges[stop]))

            old_widths[start] *= start_factor
            old_widths[stop] *= stop_factor

            f_widths = old_widths[start:stop] .* ys_old[start:stop]
            ys_new[i] = sum(f_widths)
            ys_new[i] /= sum(old_widths[start:stop])

            e_wid = old_widths[start:stop] .* σs_old[start:stop]
            σs_new[i] = sqrt(sum(e_wid.^2.0))
            σs_new[i] /= sum(old_widths[start:stop])

            old_widths[start] /= start_factor
            old_widths[stop] /= stop_factor
        end
    end
    return ys_new, σs_new
end

function convolve_gauss(xs::AA{T,1}, ys::AA{T,1}; new_res::T=1.17e5,
                        oversampling::T=1.0) where T<:AbstractFloat
    # get kernel
    # TODO: wavelength dependent kernel width???
    σ = mean(xs) / new_res / 2.354
    g(x, n) = (one(T)/(σ * sqrt(2.0 * π))) * exp(-0.5 * ((x - n)/σ)^2)
    kernel = g.(xs, xs[Int(round(length(xs)/2))])

    # pad the signal
    signal = vcat(zeros(100), ys[:,1], zeros(100))
    signal[1:100] .= first(ys[:,1])
    signal[end-100:end] .= last(ys[:,1])

    # do the convolution
    new_ys = imfilter(signal, reflect(centered(kernel./maximum(kernel))), Fill(0))
    new_ys = new_ys[101:end-100]
    new_ys ./= maximum(new_ys)

    # get wavelength grid at lower resolution
    Δlnλ = 1.0 / new_res
    lnλs = range(log(first(xs)), log(last(xs)), step=Δlnλ/oversampling)
    xs_out = exp.(lnλs)

    # flux-conserving resample the convolved data onto lower res grid
    ys_out = rebin_spectrum(xs, new_ys, xs_out)
    return xs_out, ys_out
end

function degrade_resolution(wave::AA{T,1}, spec::AA{T,1}; R_new::T=700000.0, snr::T=Inf) where T<:AF
    # get new resolution element and wavelength grid
    dλ = wave[1]/R_new
    λs = range(wave[1], wave[end], step=dλ)

    # interpolate the spectrum onto new wavelength grid
    spl = LinearInterpolation(wave, spec, extrapolation_bc=Flat())
    new_spec = spl(λs)

    # add noise and return
    return λs, add_noise(new_spec, snr)
end

function degrade_resolution(wave::AA{T,1}, spec::AA{T,2}; R_new::T=700000.0) where T<:AF
    return [degrade_resolution(wave, spec[:,i]) for i in 1:size(spec,2)]
    # return map(x -> degrade_resolution(wave, spec[:,i]), [])
end
