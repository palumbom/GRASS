function add_noise!(flux::AbstractArray{T,1}, snr::T) where T<:Real
    # TODO: snr different in continuum and line cores?
    # N = snr ^ 2.0

    # generate noise and scale to SNR
    noise = randn(length(flux))
    scale = flux ./ std(noise)
    noise .*= (scale ./ snr)
    flux .+= noise
    return nothing
end

function add_noise!(flux::AbstractArray{T,2}, snr::T) where T<:Real
    map(x -> add_noise!(x, snr), eachcol(flux))
    return nothing
end
