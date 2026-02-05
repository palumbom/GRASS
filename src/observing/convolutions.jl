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
    @assert issorted(xs_old)
    @assert issorted(xs_new)

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
    @assert issorted(xs_old)
    @assert issorted(xs_new)

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
    σ(x) = x / new_res / (2.0 * sqrt(2 * log(2)))
    g(x, n) = (one(T)/(σ(x) * sqrt(2.0 * π))) * exp(-0.5 * ((x - n)/σ(x))^2)
    kernel = g.(xs, xs[Int(round(length(xs)/2))])

    # pad the signal
    signal = vcat(zeros(100), ys[:,1], zeros(100))
    signal[1:100] .= first(ys[:,1])
    signal[101:end-100] .= ys
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

function convolve_gauss(xs::AA{T,1}, ys::AA{T,2}; new_res::T=1.17e5,
                        oversampling::T=1.0) where T<:AbstractFloat
    # get kernel
    σ(x) = x / new_res / (2.0 * sqrt(2 * log(2)))
    g(x, n) = (one(T)/(σ(x) * sqrt(2.0 * π))) * exp(-0.5 * ((x - n)/σ(x))^2)
    kernel = g.(xs, xs[Int(round(length(xs)/2))])

    # create vector to hold padded signal
    signal = zeros(200 + size(ys, 1))

    # get wavelength grid at lower resolution
    Δlnλ = 1.0 / new_res
    lnλs = range(log(first(xs)), log(last(xs)), step=Δlnλ/oversampling)
    xs_out = exp.(lnλs)

    # allocate array for output spectrum
    ys_out = zeros(length(xs_out), size(ys, 2))

    # loop over times
    for t in 1:size(ys,2)
        # fill signal vector
        signal[1:100] .= first(ys[:,t])
        signal[101:end-100] .= ys[:,t]
        signal[end-100:end] .= last(ys[:,t])

        # perform the convolution
        # TODO: fft plan (outside loop?)
        signal .= imfilter(signal, reflect(centered(kernel./maximum(kernel))), Fill(0))
        signal ./= maximum(signal[101:end-100])

        # rebin the convolved spectrum
        ys_out[:,t] .= rebin_spectrum(xs, view(signal, 101:length(signal)-100), xs_out)
    end

    return xs_out, ys_out
end

function convolve_gauss(xs::AA{T,1}, ys::AA{T,2}, kernel; new_res::T=1.17e5,
                        oversampling::T=1.0) where T<:AbstractFloat
    # create vector to hold padded signal
    signal = zeros(200 + size(ys, 1))

    # get wavelength grid at lower resolution
    Δlnλ = 1.0 / new_res
    lnλs = range(log(first(xs)), log(last(xs)), step=Δlnλ/oversampling)
    xs_out = exp.(lnλs)

    # allocate array for output spectrum
    ys_out = zeros(length(xs_out), size(ys, 2))

    # loop over times
    for t in 1:size(ys,2)
        # fill signal vector
        signal[1:100] .= first(ys[:,t])
        signal[101:end-100] .= ys[:,t]
        signal[end-100:end] .= last(ys[:,t])

        # perform the convolution
        # TODO: fft plan (outside loop?)
        signal .= imfilter(signal, reflect(centered(kernel./maximum(kernel))), Fill(0))
        signal ./= maximum(signal[101:end-100])

        # rebin the convolved spectrum
        ys_out[:,t] .= rebin_spectrum(xs, view(signal, 101:length(signal)-100), xs_out)
    end

    return xs_out, ys_out
end