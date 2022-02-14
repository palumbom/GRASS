function trim_bisector!(depth::T, wavt::AA{T,1}, bist::AA{T,1},
                        dept::AA{T,1}, widt::AA{T,1}; top::T=NaN) where T<:AF
    # choose method depending if synth line is deeper than input
    if (one(T) - depth) > minimum(bist)
        trim_bisector_chop!(depth, wavt, bist, dept, widt, top=top)
    else
        trim_bisector_scale!(depth, wavt, bist, dept, widt, top=top)
    end
    return nothing
end

function trim_bisector_chop!(depth::T, wavt::AA{T,1}, bist::AA{T,1},
                             dept::AA{T,1}, widt::AA{T,1};
                             top::T=NaN) where T<:AF
    # create interpolators
    itp1 = linear_interp(bist, wavt)
    itp2 = linear_interp(dept, widt)

    # get new grid of depths, interpolate the data, and return
    dept .= range((one(T) - depth), one(T), length=length(dept))
    wavt .= itp1.(dept)
    widt .= itp2.(dept)
    bist .= dept
    return nothing
end

function trim_bisector_scale!(depth::T, wavt::AA{T,1}, bist::AA{T,1},
                              dept::AA{T,1}, widt::AA{T,1};
                              top::T=NaN) where T<:AF
    # get new grid of depths, effectively scaling the width and bisector data
    dept .= range((one(T) - depth), one(T), length=length(dept))
    bist .= dept
    return nothing
end

