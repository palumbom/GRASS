function calc_bisector_inverse_slope(bis::AA{T,1}, int::AA{T,1}) where T<:AF
    # get total depth
    dep = one(T) - minimum(int)
    bot = minimum(int)

    # find indices
    idx10 = findfirst(x -> x .> 0.10 * dep + bot, int)
    idx40 = findfirst(x -> x .> 0.40 * dep + bot, int)
    idx55 = findfirst(x -> x .> 0.55 * dep + bot, int)
    idx90 = findfirst(x -> x .> 0.90 * dep + bot, int)

    # get v_t and v_b
    v_t = mean(bis[idx55:idx90])
    v_b = mean(bis[idx10:idx40])
    return v_t - v_b
end

function calc_bisector_inverse_slope(bis::AA{T,2}, int::AA{T,2}) where T<:AF
    out = zeros(size(bis,2))
    for i in 1:size(bis,2)
        out[i] = calc_bisector_inverse_slope(bis[:,i], int[:,i])
    end
    return out
end


function calculate_bisector_span(λrest::T, wav::AA{T,1}) where T<:AF
    minw = minimum(filter(!isnan, wav))
    return abs(minw - wav[5])/wav[5] * (c_ms)
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


function calc_bisector_uncertainty(bis::AA{T,1}, int::AA{T,1}) where T<:AF
    dF = diff(int)
    dv = diff(bis)
    return one(T)/sqrt(2.0) .* 1.0 ./ abs.(dF ./ dv)
end


function calc_line_quantity(wavs::AA{T,1}, flux::AA{T,1}; continuum::T=1.0,
                            n::Int=2, top::T=maximum(flux),
                            nflux::Int=length(flux), f::Function) where T<:Real
    # check lengths
    @assert n >= 0
    @assert length(wavs) == length(flux)

    # TODO smooth until derivative shows function is monotonic
    # perform moving average smoothing
    flux = moving_average(flux, n)
    wavs = moving_average(wavs, n)

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

    # allocate memory
    x_out = zeros(nflux)
    y_out = range(minimum(flux), top, length=nflux)

    # loop over flux values and measure quantity defined by f function
    for i in eachindex(y_out)
        lidx = findfirst(x -> x .>= y_out[i], lflux)
        ridx = findfirst(x -> x .>= y_out[i], rflux)

        if isnothing(lidx)
            lidx = length(lflux)
        end

        if isnothing(ridx)
            ridx = length(rflux)
        end

        # adjust indices to account for views
        wav_lidx = min_flux_idx - lidx
        wav_ridx = min_flux_idx + ridx

        # interpolate on left
        if lflux[lidx] != y_out[i]
            w2 = (y_out[i] - lflux[lidx-1]) / (lflux[lidx] - lflux[lidx-1])
            w1 = 1.0 - w2
            lwav = lwavs[lidx-1] * w1 + lwavs[lidx] * w2
        else
            lwav = lwavs[ridx]
        end

        # interpolate on right
        if rflux[ridx] != y_out[i]
            w2 = (y_out[i] - rflux[ridx-1]) / (rflux[ridx] - rflux[ridx-1])
            w1 = 1.0 - w2
            rwav = rwavs[ridx-1] * w1 + rwavs[ridx] * w2
        else
            rwav = rwavs[ridx]
        end

        # calculate and assign quantity value
        x_out[i] = f(lwav, rwav)
    end
    return x_out, y_out
end


function calc_width_function(wavs::AA{T,1}, flux::AA{T,1}; kwargs...) where T<:Real
    # check lengths
    @assert length(wavs) == length(flux)

    # define the function to calculate widths
    f = (x, y) -> (y - x)
    wid, dep = calc_line_quantity(wavs, flux, f=f; kwargs...)
    return dep, wid
end


function calc_width_function(wavs::AA{T,2}, flux::AA{T,2}; kwargs...) where T<:Real
    f = (x,y) -> calc_width_function(x, y; kwargs...)
    out = map(f, eachcol(wavs), eachcol(flux))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end


function calc_bisector(wavs::AA{T,1}, flux::AA{T,1}; kwargs...) where T<:Real
    # check lengths
    @assert length(wavs) == length(flux)

    # define the function to calculate widths
    f = (x, y) -> (y + x) / 2.0
    wav, bis = calc_line_quantity(wavs, flux, f=f; kwargs...)
    return wav, bis
end


function calc_bisector(wavs::AA{T,2}, flux::AA{T,2}; kwargs...) where T<:Real
    f = (x,y) -> calc_bisector(x, y; kwargs...)
    out = map(f, eachcol(wavs), eachcol(flux))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end

function calc_bisector(wavs::AA{T,1}, flux::AA{T,2}; kwargs...) where T<:Real
    f = y -> calc_bisector(wavs, y; kwargs...)
    out = map(f, eachcol(flux))
    return cat([x[1] for x in out]..., dims=2), cat([x[2] for x in out]..., dims=2)
end

function calc_bisector_cubic(wavs::AA{T,1}, flux::AA{T,1}) where T<:Real
    # find index of minimum
    idx_min = argmin(flux)

    # take difference quotient
    dfdλ = diff(flux)./diff(wavs)

    # starting from middle, find last element where deriv. is pos/neg.
    idx1 = idx_min - findlast(view(dfdλ, idx_min:-1:1) .< 0.0) + 1
    idx2 = findlast(view(dfdλ, idx_min:length(dfdλ)) .> 0.0) + idx_min - 1

    # construct the interpolants
    itp1 = cubic_interp(view(flux, idx_min:-1:idx1), view(wavs, idx_min:-1:idx1))
    itp2 = cubic_interp(view(flux, idx_min:idx2), view(wavs, idx_min:idx2))

    # get the bisector
    flux_out = range(minimum(flux), maximum(flux) - 0.01, length=50)
    bis_out = (itp2.(flux_out) .+ itp1.(flux_out))./2

    return bis_out, flux_out
end
