function linear_interp_gpu(xs, ys, bc)


    return
end

function trim_bisector_chop_gpu!(depth, wavt, bist, dept, widt, top)
    # replace spurious measurements at top of bisector
    ind1 = searchsortedfirst(bist, one(T) - depth)
    if !isnan(top)
        ind2 = searchsortedfirst(bist, top)
        wavt[ind2:end] .= wavt[ind2]
    end

    # set up interpolant
    itp1 = linear_interp(view(bist, ind1:length(bist)), view(wavt, ind1:length(wavt)))

    # get new grid of depths, interpolate the data, and return
    dept .= range((one(T) - depth), one(T), length=length(dept))
    wavt .= itp1.(dept)
    bist .= dept
    return nothing
end


function line_loop_gpu(prof, mid, depth, rot_shift, conv_blueshift, lambdas, wsp, top)


    return
end

xs = CUDA.range(-5, 5, length=100)
ys = xs.^2

function linear_interp_gpu(out, new_xs, xs, ys)
    # get GPU dims
    i = threadIdx().x + blockDim().x * (blockIdx().x-1)
    # j = threadIdx().y + blockDim().y * (blockIdx().y-1)

    # get length of data
    n = length(xs)

    if i < n
        j = CUDA.searchsortedfirst(xs, new_xs[i]) - 1
        j0 = CUDA.clamp(j, CUDA.firstindex(ys), CUDA.lastindex(ys))
        j1 = CUDA.clamp(j+1, CUDA.firstindex(ys), CUDA.lastindex(ys))

        out[i] = ys[j0] + (ys[j1] - ys[j0]) * (new_xs[i] - xs[j0]) / (xs[j1] - xs[j0])
    end
    return
end

xs = CuArray(range(-5, 5, length=100))
ys = CuArray()
