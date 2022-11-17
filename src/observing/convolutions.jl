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
function rebin_spectrum(xs_old, ys_old, xs_new)
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

function convolve_gauss(xs::AA{T,1}, ys::AA{T,1}; new_res::T=1.17e5,
                        oversampling::T=1.0) where T<:AbstractFloat
    # TODO: input spectrum is not infinitely sharp
    # new res is actually convolution res
    σ = (5.0/7.0) * 5434.5232 / new_res / 2.354
    g(x, n) = (one(T)/(σ * sqrt(2.0 * π))) * exp(-0.5 * ((x - n)/σ)^2)

    # pad x array to deal with edges
    # TODO: make padding smarter
    xstep = minimum(diff(xs))
    newxs = range(first(xs)-100*xstep, last(xs)+100*xstep, step=xstep)
    newxs = vcat(range(first(xs)-100*xstep, first(xs)-xstep, step=xstep),
                 xs,
                 range(last(xs)+xstep, last(xs)+100*xstep, step=xstep))

    # find matching indices
    ind1 = searchsortedfirst(newxs, xs[1])
    ind2 = searchsortedfirst(newxs, xs[end])

    # pad y array
    newys = similar(newxs)
    newys[1:ind1-1] .= first(ys)
    newys[ind1:ind2] .= ys
    newys[ind2+1:end] .= last(ys)

    # perform convolution
    conv = similar(newxs)
    eval_g = zeros(size(newys))
    for i in eachindex(conv)
        eval_g .= g.(newxs, newxs[i])
        conv[i] = sum(newys .* eval_g) / sum(eval_g)
    end

    # get wavelength grid at lower resolution
    Δlnλ = 1.0 / (new_res * oversampling)
    lnλs = range(log(first(xs)), log(last(xs)), step=Δlnλ)
    xs_out = exp.(lnλs)

    # re-sample convolved data onto lower res grid (preserving flux)
    ys_out = rebin_spectrum(newxs, newys, xs_out)
    return xs_out, ys_out
end

function convolve_gauss(xs::AA{T,1}, ys::AA{T,2}; new_res::T=1.17e5) where T<:AbstractFloat
    conv = zeros(length(xs), size(ys,2))
    for i in 1:size(ys,2)
        conv[:, i] = convolve_gauss(xs, ys[:,i], new_res=new_res)
    end
    return conv
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
