function add_noise(spec::Array{T,1}, snr::T) where T<:Real
    # generate noise and scale to SNR
    noise = randn(length(spec))
    scale = spec ./ std(noise)
    noise .*= (scale ./ snr)
    spec .+= noise
    return spec
end

function add_noise(spec::Array{T,2}, snr::T) where T<:Real
    # loop over each time slice adding noise
    spec_noise = zeros(size(spec))
    for t in 1:size(spec,2)
        spec_noise[:,t] = add_noise(spec[:,t], snr)
    end
    return spec_noise
end
