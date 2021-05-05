function convolve_gauss(xs::AA{T,1}, ys::AA{T,1}; new_res::T=1.17e5) where T<:AbstractFloat
    # TODO: input spectrum is not infinitely sharp
    # new res is actually convolution res
    σ = (5.0/7.0) * 5434.5232 / new_res / 2.354
    g(x, n) = (one(T)/(σ * sqrt(2.0 * π))) * exp(-0.5 * ((x - n)/σ)^2)

    # pad x array to deal with edges
    # TODO: make padding smarter
    xstep = step(xs)
    newxs = range(xs[1]-100*xstep, xs[end]+100*xstep, step=xstep)

    # find matching indices
    ind1 = searchsortedfirst(newxs, xs[1])
    ind2 = searchsortedfirst(newxs, xs[end])

    # pad y array
    newys = similar(newxs)
    newys[1:ind1-1] .= 1.0
    newys[ind1:ind2] .= ys
    newys[ind2+1:end] .= 1.0

    # perform convolution
    # TODO: do convolution on different grid than newxs
    conv = similar(newxs)
    for i in eachindex(conv)
        conv[i] = sum(newys .* g.(newxs, newxs[i]))
    end

    # normalize and return sliced convolution
    conv ./= maximum(conv)
    return conv[ind1:ind2]
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
