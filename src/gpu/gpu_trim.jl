function trim_bisector_chop_gpu(depth, wavall_out, bisall_out, depall_out, widall_out,
                                wavall_in, bisall_in, depall_in, widall_in, top)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # loop over epochs of bisectors
    for i in idx:sdx:CUDA.size(wavall_in,2)
        # get views of the correct time slice in input
        wavt = CUDA.view(wavall_in, :, i)
        bist = CUDA.view(bisall_in, :, i)
        dept = CUDA.view(depall_in, :, i)
        widt = CUDA.view(widall_in, :, i)

        # get views of the correct time slice in out
        wavt_out = CUDA.view(wavall_out, :, i)
        bist_out = CUDA.view(bisall_out, :, i)
        dept_out = CUDA.view(depall_out, :, i)
        widt_out = CUDA.view(widall_out, :, i)

        # find first index above specified depth
        ind1 = CUDA.searchsortedfirst(bist, 1.0 - depth)

        # TODO this will kill the code (I think)
        if !CUDA.isnan(top)
            ind2 = CUDA.searchsortedfirst(bist, top)
            wavt[ind2:end] .= wavt[ind2]
        end

        # get knots
        xs = CUDA.view(bist, ind1:CUDA.length(bist))
        ys = CUDA.view(wavt, ind1:CUDA.length(wavt))

        # TODO: writing to same memory other threads are reading which is weird
        # xs = bist
        # ys = dept

        # loop over the length of the bisector
        step = depth/(CUDA.length(dept) - 1)
        CUDA.@assert(typeof(step) == Float64)
        for j in idy:sdy:CUDA.size(wavall_in,1)
            # set the new depth value
            new_dept = (1.0 - depth) + (j-1) * step

            # interpolate to get the new wavelength value
            if (((new_dept < CUDA.first(xs)) || (new_dept > CUDA.last(xs))) && !CUDA.isnan(top))
                wavt_out[j] = top
            elseif new_dept <= CUDA.first(xs)
                wavt_out[j] = CUDA.first(ys)
            elseif new_dept >= CUDA.last(xs)
                wavt_out[j] = CUDA.last(ys)
            else
                k = CUDA.searchsortedfirst(xs, new_dept) - 1
                k0 = CUDA.clamp(k, CUDA.firstindex(ys), CUDA.lastindex(ys))
                k1 = CUDA.clamp(k+1, CUDA.firstindex(ys), CUDA.lastindex(ys))
                wavt_out[j] = ys[k0] + (ys[k1] - ys[k0]) * (new_dept - xs[k0]) / (xs[k1] - xs[k0])
            end

            # assign bisector fluxes from dept
            bist_out[j] = new_dept
            dept_out[j] = new_dept # TODO this is unnecessary
        end
    end
    return nothing
end
