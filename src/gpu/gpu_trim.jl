# TODO fix this function for the case where template is shallower than synthetic
function trim_bisector_gpu(depth, wavall_out, depall_out, wavall_in, depall_in)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # loop over epochs of bisectors
    for i in idx:sdx:CUDA.size(wavall_in, 2)
        # get views of the correct time slice in input
        wavt_in = CUDA.view(wavall_in, :, i)
        dept_in = CUDA.view(depall_in, :, i)

        # get views of the correct time slice in output
        wavt_out = CUDA.view(wavall_out, :, i)
        dept_out = CUDA.view(depall_out, :, i)

        # loop over the length of the bisector
        step = depth/(CUDA.length(dept_in) - 1)
        for j in idy:sdy:CUDA.size(wavall_in, 1)
            # set the new depth value
            new_dept = (1.0 - depth) + (j-1) * step

            # if depth is greater than input depth, stretch the bisector
            if (1.0 - depth) < CUDA.first(dept_in)
                @inbounds dept_out[j] = new_dept
            else
                # interpolate to get the new wavelength value
                if new_dept <= CUDA.first(dept_in)
                    @inbounds wavt_out[j] = CUDA.first(wavt_in)
                elseif new_dept >= CUDA.last(dept_in)
                    @inbounds wavt_out[j] = CUDA.last(wavt_in)
                else
                    k = CUDA.searchsortedfirst(dept_in, new_dept) - 1
                    k0 = CUDA.clamp(k, CUDA.firstindex(wavt_in), CUDA.lastindex(wavt_in))
                    k1 = CUDA.clamp(k+1, CUDA.firstindex(wavt_in), CUDA.lastindex(wavt_in))
                    @inbounds wavt_out[j] = (wavt_in[k0] * (dept_in[k1] - new_dept) + wavt_in[k1] * (new_dept - dept_in[k0])) / (dept_in[k1] - dept_in[k0])
                end

                # assign bisector fluxes from dept
                @inbounds dept_out[j] = new_dept
            end
        end
    end
    return nothing
end
