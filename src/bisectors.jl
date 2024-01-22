"""
From Queloz et al. 2001
"""
function calc_bisector_inverse_slope(bis::AA{T,1}, int::AA{T,1}; b1::T=0.10,
                                     b2::T=0.4, b3::T=0.55, b4::T=0.90) where T<:AF
    @assert maximum(int) <= 1.0
    @assert minimum(int) >= 0.0

    # get total depth
    dep = one(T) - minimum(int)
    bot = minimum(int)
    top = one(T)

    # find indices
    idx10 = findfirst(x -> x .> top - b1 * dep, int)
    idx40 = findfirst(x -> x .> top - b2 * dep, int)
    idx55 = findfirst(x -> x .> top - b3 * dep, int)
    idx90 = findfirst(x -> x .> top - b4 * dep, int)

    # get v_t and v_b
    v_t = mean(view(bis, idx40:idx10))
    v_b = mean(view(bis, idx90:idx55))
    return v_t - v_b
end

function calc_bisector_inverse_slope(bis::AA{T,2}, int::AA{T,2}; b1::T=0.10,
                                     b2::T=0.4, b3::T=0.55, b4::T=0.90) where T<:AF
    out = zeros(size(bis,2))
    for i in 1:size(bis,2)
        out[i] = calc_bisector_inverse_slope(bis[:,i], int[:,i], b1=b1, b2=b2, b3=b3, b4=b4)
    end
    return out
end

"""
From various (TODO)
"""
function calc_bisector_span(bis::AA{T,1}, int::AA{T,1}) where T<:AF
    @assert maximum(int) <= 1.0
    @assert minimum(int) >= 0.0

    blue = minimum(bis[3:end])
    core = mean(bis[3:5])
    return core - blue
end


function calc_bisector_span(bis::AA{T,2}, int::AA{T,2}) where T<:AF
    out = zeros(size(bis,2))
    for i in 1:size(bis,2)
        out[i] = calc_bisector_span(bis[:,i], int[:,i])
    end
    return out
end

"""
From Dall et al. 2006
"""
function calc_bisector_slope(bis::AA{T,1}, int::AA{T,1}) where T<:AF
    @assert maximum(int) <= 1.0
    @assert minimum(int) >= 0.0

    # get total depth
    dep = one(T) - minimum(int)
    bot = minimum(int)
    top = one(T)

    # find indices
    idx25 = findfirst(x -> x .> top - 0.25 * dep, int)
    idx80 = findfirst(x -> x .> top - 0.80 * dep, int)

    # get views
    bis_view = view(bis, idx80:idx25)
    int_view = view(int, idx80:idx25)

    # linear fit
    pfit = Polynomials.fit(bis_view, int_view, 1)

    # return inverse slope
    return -coeffs(pfit)[2]
end

function calc_bisector_slope(bis::AA{T,2}, int::AA{T,2}) where T<:AF
    out = zeros(size(bis,2))
    for i in 1:size(bis,2)
        out[i] = calc_bisector_slope(bis[:,i], int[:,i])
    end
    return out
end

"""
From Dall et al. 2006
"""
function calc_bisector_bottom(bis::AA{T,1}, int::AA{T,1}, rv::T) where T<:AF
    @assert maximum(int) <= 1.0
    @assert minimum(int) >= 0.0

    # return inverse slope
    return mean(view(bis, 2:5)) - rv
end

function calc_bisector_bottom(bis::AA{T,2}, int::AA{T,2}, rvs::AA{T,1}) where T<:AF
    out = zeros(size(bis,2))
    for i in 1:size(bis,2)
        out[i] = calc_bisector_bottom(bis[:,i], int[:,i], rvs[i])
    end
    return out
end


"""
From Dall et al. 2006
"""
function calc_bisector_curvature(bis::AA{T,1}, int::AA{T,1}) where T<:AF
    @assert maximum(int) <= 1.0
    @assert minimum(int) >= 0.0

    # get total depth
    dep = one(T) - minimum(int)
    bot = minimum(int)
    top = one(T)

    # find indices
    idx20 = findfirst(x -> x .> top - 0.20 * dep, int)
    idx30 = findfirst(x -> x .> top - 0.30 * dep, int)
    idx40 = findfirst(x -> x .> top - 0.40 * dep, int)
    idx55 = findfirst(x -> x .> top - 0.55 * dep, int)
    idx75 = findfirst(x -> x .> top - 0.75 * dep, int)
    idx100 = 1

    # take views
    v3 = mean(view(bis, idx100:idx75))
    v2 = mean(view(bis, idx55:idx40))
    v1 = mean(view(bis, idx30:idx20))

    return (v3 - v2) - (v2 - v1)
end

function calc_bisector_curvature(bis::AA{T,2}, int::AA{T,2}) where T<:AF
    out = zeros(size(bis,2))
    for i in 1:size(bis,2)
        out[i] = calc_bisector_curvature(bis[:,i], int[:,i])
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

    # find next-to-minimum flux value
    idx_sort = sortperm(flux)
    min_meas = flux[idx_sort[5]]

    # allocate memory
    x_out = zeros(nflux)
    y_out = range(min_meas, top, length=nflux)

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

    # define the function to calculate bisectors
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

function calc_bisector_cubic(wavs::AA{T,1}, flux::AA{T,1}; top::T=maximum(flux),
                             nflux::Int=length(flux)) where T<:Real
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
    flux_out = range(minimum(flux), top, length=nflux)
    bis_out = (itp2.(flux_out) .+ itp1.(flux_out))./2

    return bis_out, flux_out
end
