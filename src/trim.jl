function trim_bisector_chop!(depth::T, wavt::AA{T,1}, bist::AA{T,1},
                             dept::AA{T,1}, widt::AA{T,1};
                             top::T=NaN) where T<:AF
    # replace spurious measurements at top of bisector
    ind1 = searchsortedfirst(bist, one(T) - depth)
    if !isnan(top)
        ind2 = searchsortedfirst(bist, top)
        wavt[ind2:end] .= wavt[ind2]
    end

    # get new grid of depths
    dept .= range((one(T) - depth), one(T), length=length(dept))
    bist .= dept

    # spline interpolate onto same-length grid
    itp1 = extrapolate(interpolate!(wavt, BSpline(Linear())), Flat())
    wavt .= itp1.(range(ind1, length(dept), length=length(dept)))
    return nothing
end
