function trim_bisector_chop_gpu!(depth, wavall, bisall, depall, widall, top)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # loop over epochs of bisectors
    for i in idx:sdx:CUDA.size(wavall,2)
        # get views of the correct slice
        wavt = CUDA.view(wavall, :, i)
        bist = CUDA.view(bisall, :, i)
        dept = CUDA.view(depall, :, i)
        widt = CUDA.view(widall, :, i)

        # find first index above specified depth
        ind1 = CUDA.searchsortedfirst(bist, 1.0 - depth)

        # TODO this will kill the code
        if !CUDA.isnan(top)
            ind2 = CUDA.searchsortedfirst(bist, top)
            wavt[ind2:end] .= wavt[ind2]
        end

        # get knots
        xs = CUDA.view(bist, ind1:CUDA.length(bist))
        ys = CUDA.view(wavt, ind1:CUDA.length(wavt))

        # loop over the length of the bisector
        step = depth/(CUDA.length(dept) - 1)
        for j in idy:sdy:CUDA.size(wavall,1)
            # set the new depth value=
            new_dept = (1.0 - depth) + (j-1) * step

            # interpolate to get the new wavelength value
            if (((new_dept < CUDA.first(xs)) | (new_dept > CUDA.last(xs))) & !CUDA.isnan(top))
                wavt[j] = top
            elseif new_dept <= CUDA.first(xs)
                wavt[j] = CUDA.first(ys)
            elseif new_dept >= CUDA.last(xs)
                wavt[j] = CUDA.last(ys)
            else
                k = CUDA.searchsortedfirst(xs, new_dept) - 1
                k0 = CUDA.clamp(k, CUDA.firstindex(ys), CUDA.lastindex(ys))
                k1 = CUDA.clamp(k+1, CUDA.firstindex(ys), CUDA.lastindex(ys))
                wavt[j] = ys[k0] + (ys[k1] - ys[k0]) * (new_dept - xs[k0]) / (xs[k1] - xs[k0])
            end

            # assign bisector fluxes from dept
            bist[j] = new_dept
            dept[j] = new_dept
        end
    end
    return nothing
end
