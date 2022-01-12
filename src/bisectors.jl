function calculate_bisector_span(λrest::T, wav::AA{T,1}) where T<:AF
    minw = minimum(filter(!isnan, wav))
    return abs(minw - wav[1])/wav[1] * (c /100.0)
end

function calculate_bisector_span(λrest::T, wav::AA{T,2}) where T<:AF
    out = zeros(size(wav,2))
    for i in eachindex(out)
        out[i] = calculate_bisector_span(λrest, wav[:,i])
    end
    return out
end

function calculate_bisector_extreme(λrest::T, wav::AA{T,2}, bis::AA{T,2}) where T<:AF
    out = zeros(size(wav,2))
    for i in eachindex(out)
        out[i] = calculate_bisector_extreme(λrest, wav[:,i], bis[:,i])
    end
    return out
end

function calculate_bisector_extreme(λrest::T, wav::AA{T,1}, bis::AA{T,1}) where T<:AF
    ind1 = searchsortednearest(wav, 0.4)
    ind2 = searchsortednearest(wav, 0.8)
    wpos = mean(wav[ind1:ind2])
    return (λrest - wpos)/λrest * c/100.0
end

function calculate_bisector_lslope(λrest::T, wav::AA{T,1}, bis::AA{T,1}) where T<:AF
    dλ = (minimum(wav[.!isnan.(wav)]) - wav[1])/λrest * c/100.0
    ind1 = findfirst(wav .== minimum(wav[.!isnan.(wav)]))

    dF = bis[ind1] - bis[1]
    return dF/dλ
end

function calculate_bisector_lslope(λrest::T, wav::AA{T,2}, bis::AA{T,2}) where T<:AF
    out = zeros(size(wav,2))
    for i in eachindex(out)
        out[i] = calculate_bisector_lslope(λrest, wav[:,i], bis[:,i])
    end
    return out
end

function measure_bisector(xs::AA{T,1}, ys::AA{T,1}; interpolate::Bool=true,
                          top::T=0.99, len::Integer=100) where T<:AF
    if interpolate
        return measure_bisector_interpolate(xs, ys, top=top, len=len)
    else
        return measure_bisector_loop(xs, ys, top=top, len=len)
    end
end

function measure_bisector_interpolate(xs::AA{T,1}, ys::AA{T,1}; top::T=0.99,
                                      len::Integer=100, max_loop::Int=20) where T<:AF
    # check lengths and normalization
    @assert length(xs) == length(ys)

    # normalize the spec, find bottom of line
    ys ./= maximum(ys)
    botind = argmin(ys)
    depths = range(ys[botind], top, length=len)

    # find left and right halves
    lind = findfirst(ys .< top)
    rind = findlast(ys .< top)
    lspec = reverse(ys[lind:botind])
    rspec = ys[botind:rind]
    lwav = reverse(xs[lind:botind])
    rwav = xs[botind:rind]

    # make sure lspec is sorted
    num_loop = 0
    unsorted = !issorted(lspec)
    while unsorted
        num_loop += 1
        ldiff = diff(lspec)
        inds = BitArray(vcat(0, ldiff .< 0))
        lspec[inds] .= NaN
        itp = LinearInterpolation(reverse(lwav[.!inds]),
                                  reverse(lspec[.!inds]),
                                  extrapolation_bc=Flat())
        lspec .= itp.(lwav)
        unsorted = !issorted(lspec)
        if num_loop > max_loop
            break
        end
    end

    # make sure rspec is sorted
    num_loop = 0
    unsorted = !issorted(rspec)
    while unsorted
        num_loop += 1
        rdiff = diff(rspec)
        inds = BitArray(vcat(0, rdiff .< 0))
        rspec[inds] .= NaN
        itp = LinearInterpolation(rwav[.!inds], rspec[.!inds],
                                  extrapolation_bc=Flat())
        rspec .= itp.(rwav)
        unsorted = !issorted(rspec)
        if num_loop > max_loop
            break
        end
    end

    # interpolate wavelengths onto intensity grid
    lspline = LinearInterpolation(lspec, lwav, extrapolation_bc=Flat())
    rspline = LinearInterpolation(rspec, rwav, extrapolation_bc=Flat())
    wavs = (lspline(depths) .+ rspline(depths)) ./ 2.0
    return wavs, depths
end

function measure_bisector_loop(xs::AA{T,1}, ys::AA{T,1}; top::T=0.99,
                               len::Integer=100) where T<:AF
    # normalize the spec, find bottom of line
    ys ./= maximum(ys)

    # assign depths to measure bisector at
    dep = range(one(T)-minimum(ys)-0.01, one(T) - top, length=len)

    # set iterators
    nccf = Int(length(xs) ÷ 2)
    L = nccf
    R = nccf

    # allocate memory
    xL = zeros(len)
    xR = zeros(len)
    wav = zeros(len)

    # loop over depths
    for d in eachindex(dep)
        y = one(T) - dep[d]
        while((ys[L] < y) & (L > 0))
            L -= 1
        end

        while ((ys[R] < y) & (R < length(xs)))
            R += 1
        end

        if ((y > maximum(ys[1:nccf])) | (y > maximum(ys[nccf+1:end])))
            L = 0
            R = length(xs)
        end

        if L == 0
            xL[d] = xL[d-1]
        else
            mL = (xs[L+1] - xs[L]) / (ys[L+1] - ys[L])
            xL[d] = xs[L] + mL * (y - ys[L])
        end

        if R == length(xs)
            xR[d] = xR[d-1]
        else
            mR = (xs[R-1] - xs[R]) / (ys[R-1] - ys[R])
            xR[d] = xs[R] + mR * (y - ys[R])
        end
        wav[d] = (xL[d] + xR[d]) / 2.0
    end
    return wav, one(T) .- dep
end

function bisector_uncertainty(wav::AA{T,1}, bis::AA{T,1}) where T<:AF
    dF = diff(bis)
    dv = diff(wav)
    return one(T)/sqrt(2.0) .* dF ./ abs.(dF ./ dv)
end


# TODO preliminary
function calc_bisector(wavs::AbstractArray{T,1}, flux::AbstractArray{T,1}; continuum::T=1.0) where T<:Real
    # check lengths
    @assert length(wavs) == length(flux)

    # get min and max flux, depth of line
    min_flux_idx = argmin(flux)
    min_flux = minimum(flux)
    max_flux = continuum
    depth = 1.0 - min_flux/max_flux # TODO check

    # get views on either side of line
    lflux = view(flux, min_flux_idx:-1:1)
    rflux = view(flux, min_flux_idx:length(flux))
    lwavs = view(wavs, min_flux_idx:-1:1)
    rwavs = view(wavs, min_flux_idx:length(wavs))

    # range of fluxes to measure bisector at
    bis_flux = range(minimum(flux), maximum(flux), length=100)

    # allocate memory for wavelength values
    bis_wavs = similar(bis_flux)

    # loop over flux values and measure bisector
    for i in eachindex(bis_flux)
        lidx = searchsortedfirst(lflux, bis_flux[i])
        ridx = searchsortedfirst(rflux, bis_flux[i])
        if (lidx > length(bis_flux)) || (ridx > length(bis_flux))
            bis_wavs[i:end] .= NaN
            break
        else
            # adjust indices to account for views
            wav_lidx = min_flux_idx - lidx
            wav_ridx = min_flux_idx + ridx

            # interpolate on left
            if lflux[lidx] != bis_flux[i]
                w2 = (bis_flux[i] - lflux[lidx-1]) / (lflux[lidx] - lflux[lidx-1])
                w1 = 1.0 - w2
                lwav = lwavs[lidx-1] * w1 + lwavs[lidx] * w2
            else
                lwav = lwavs[ridx]
            end

            # interpolate on right
            if rflux[ridx] != bis_flux[i]
                w2 = (bis_flux[i] - rflux[ridx-1]) / (rflux[ridx] - rflux[ridx-1])
                w1 = 1.0 - w2
                rwav = rwavs[ridx-1] * w1 + rwavs[ridx] * w2
            else
                rwav = rwavs[ridx]
            end

            # bis_wavs[i] = (wavs[ridx] + wavs[lidx]) ./ 2.0
            bis_wavs[i] = (rwav + lwav) ./ 2.0
        end
    end
    return bis_wavs, bis_flux
end
