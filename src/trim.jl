function trim_bisector_chop!(depth::T, wavt::AA{T,1}, bist::AA{T,1},
                             dept::AA{T,1}, widt::AA{T,1};
                             top::T=NaN) where T<:AF
    # create interpolator
    itp1 = linear_interp(bist, wavt)

    # get new grid of depths, interpolate the data, and return
    dept .= range((one(T) - depth), one(T), length=length(dept))
    wavt .= itp1.(dept)
    bist .= dept
    return nothing
end
