"""
function linear_interp_gpu(out, new_xs, xs, ys, bc)
    # perform the interpolation
    n = CUDA.length(new_xs)

    for i in 1:CUDA.length(new_xs)
        if (((new_xs[i] < CUDA.first(xs)) | (new_xs[i] > CUDA.last(xs))) & !CUDA.isnan(bc))
            out[i] = bc
        elseif new_xs[i] <= CUDA.first(xs)
            out[i] = CUDA.first(ys)
        elseif new_xs[i] >= CUDA.last(xs)
            out[i] = CUDA.last(ys)
        else
            j = CUDA.searchsortedfirst(xs, new_xs[i]) - 1
            j0 = CUDA.clamp(j, CUDA.firstindex(ys), CUDA.lastindex(ys))
            j1 = CUDA.clamp(j+1, CUDA.firstindex(ys), CUDA.lastindex(ys))
            out[i] = ys[j0] + (ys[j1] - ys[j0]) * (new_xs[i] - xs[j0]) / (xs[j1] - xs[j0])
        end
    end
    return nothing
end

function linear_interp_mult_gpu(newx, xs, ys, bc)
    if (((newx < CUDA.first(xs)) | (newx > CUDA.last(xs))) & !CUDA.isnan(bc))
        return bc
    elseif newx <= CUDA.first(xs)
        return CUDA.first(ys)
    elseif newx >= CUDA.last(xs)
        return CUDA.last(ys)
    else
        m = CUDA.searchsortedfirst(xs, newx) - 1
        m0 = CUDA.clamp(m, CUDA.firstindex(ys), CUDA.lastindex(ys))
        m1 = CUDA.clamp(m+1, CUDA.firstindex(ys), CUDA.lastindex(ys))
        return (ys[m0] + (ys[m1] - ys[m0]) * (newx - xs[m0]) / (xs[m1] - xs[m0]))
    end
end
"""
