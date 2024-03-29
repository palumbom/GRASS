function trim_bisector_gpu!(depth, variability, depcontrast, lenall, bisall_out,
                            intall_out, widall_out, bisall_in, intall_in, widall_in)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y
    idz = threadIdx().z + blockDim().z * (blockIdx().z-1)
    sdz = blockDim().z * gridDim().z

    # loop over disk positions for bisectors
    for i in idx:sdx:CUDA.length(lenall)
        # get depth to trim to
        dtrim = depth * depcontrast[i]

        # loop over epochs of bisectors
        for j in idy:sdy:CUDA.size(bisall_in, 2)
            # don't bother trimming nothing
            if j > lenall[i]
                continue
            end

            # get views of the correct time slice in input
            bist_in = CUDA.view(bisall_in, :, j, i)
            intt_in = CUDA.view(intall_in, :, j, i)

            # get views of the correct time slice in output
            bist_out = CUDA.view(bisall_out, :, j, i)
            intt_out = CUDA.view(intall_out, :, j, i)

            # get step size for loop over length of bisector
            step = dtrim/(CUDA.length(intt_in) - 1)

            if variability
                # set up interpolator
                itp = linear_interp_gpu(intt_in, bist_in)

                # loop over the length of the bisector
                for k in idz:sdz:CUDA.size(bisall_in, 1)
                    new_intt = (1.0 - dtrim) + (k-1) * step
                    if (1.0 - dtrim) >= CUDA.first(intt_in)
                        @inbounds bist_out[k] = itp(new_intt)
                    end
                    @inbounds intt_out[k] = new_intt
                end
            else
                for k in idz:sdz:CUDA.size(bisall_in, 1)
                    @inbounds intt_out[k] = (1.0 - dtrim) + (k-1) * step
                    @inbounds bist_out[k] = 0.0
                    @inbounds widall_out[k,j,i] = widall_in[k,1,i]
                end
            end
        end
    end
    return nothing
end
