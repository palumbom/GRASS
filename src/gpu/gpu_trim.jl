function trim_bisector_chop_gpu!(depth, wavt, bist, dept, widt, top)
    # replace spurious measurements at top of bisector
    ind1 = CUDA.searchsortedfirst(bist, 1.0 - depth)

    # TODO this will kill the code
    if !CUDA.isnan(top)
        ind2 = CUDA.searchsortedfirst(bist, top)
        wavt[ind2:end] .= wavt[ind2]
    end

    # get knots
    xs = CUDA.view(bist, ind1:CUDA.length(bist))
    ys = CUDA.view(wavt, ind1:CUDA.length(wavt))

    # get new grid of depths
    step = depth/(CUDA.length(dept) - 1)
    for i in 1:CUDA.length(dept)
        dept[i] = (1.0 - depth) + (i-1) * step
    end

    # do the interpolation
    linear_interp_gpu(wavt, dept, xs, ys, NaN)

    # now assign bisector fluxes from dept
    for i in 1:CUDA.length(bist)
        bist[i] = dept[i]
    end
    return nothing
end
