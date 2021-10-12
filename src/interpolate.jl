function linear_interp(xs::AbstractArray{T,1}, ys::AbstractArray{T,1}; bc::T=1.0) where T<:Float64
    function f(x)
        # if ((x <= first(xs)) | (x >= last(xs)))
        #     return bc
        if x <= first(xs)
            return first(ys)
        elseif x >= last(xs)
            return last(ys)
        else
            i = searchsortedfirst(xs, x) - 1
            i0 = clamp(i, firstindex(ys), lastindex(ys))
            i1 = clamp(i+1, firstindex(ys), lastindex(ys))
            return ys[i0] + (ys[i1] - ys[i0]) * (x - xs[i0]) / (xs[i1] - xs[i0])
        end
    end
    return f
end
