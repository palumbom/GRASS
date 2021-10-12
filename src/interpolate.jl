function linear_interp(xs::AbstractArray{T,1}, ys::AbstractArray{T,1}; bc::T=NaN) where T<:Float64
    function f(x)
        if (((x <= first(xs)) | (x >= last(xs))) & !isnan(bc))
            return bc
        elseif x <= first(xs)
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
