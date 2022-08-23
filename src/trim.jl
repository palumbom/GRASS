function trim_bisector!(depth::T, bist::AA{T,1}, intt::AA{T,1}) where T<:AF
    # choose method depending if synth line is deeper than input
    if (one(T) - depth) > minimum(bist)
        trim_bisector_chop!(depth, bist, intt)
    else
        trim_bisector_scale!(depth, bist, intt)
    end
    return nothing
end

function trim_bisector_chop!(depth::T, bist::AA{T,1}, intt::AA{T,1}) where T<:AF
    # create interpolators
    itp1 = linear_interp(intt, bist)

    # get new grid of depths, interpolate the data, and return
    # TODO figure out memory with new_int
    new_int = range((one(T) - depth), one(T), length=length(intt))
    bist .= itp1.(new_int)
    intt .= new_int
    return nothing
end

function trim_bisector_scale!(depth::T, bist::AA{T,1}, intt::AA{T,1}) where T<:AF
    # get new grid of depths, effectively scaling the width and bisector data
    intt .= range((one(T) - depth), one(T), length=length(intt))
    return nothing
end

