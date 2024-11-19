function trim_bisector_gpu!(depth, variability, depcontrast, lenall, bisall_out,
                            intall_out, widall_out, bisall_in, intall_in, widall_in,
                            bisall_mean, intall_mean, widall_mean)
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

            # get views of the correct time slice in output
            bist_out = CUDA.view(bisall_out, :, j, i)
            intt_out = CUDA.view(intall_out, :, j, i)
            widt_out = CUDA.view(widall_out, :, j, i)

            if variability
                # get views of the correct time slice in input
                bist_in = CUDA.view(bisall_in, :, j, i)
                intt_in = CUDA.view(intall_in, :, j, i)

                # get step size for loop over length of bisector
                step = dtrim/(CUDA.length(intt_in) - 1)

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
                # get views of the correct time slice in input
                bist_in = CUDA.view(bisall_mean, :, i)
                intt_in = CUDA.view(intall_mean, :, i)
                widt_in = CUDA.view(widall_mean, :, i)
                
                # get step size for loop over length of bisector
                step = dtrim/(CUDA.length(intt_in) - 1)

                # set up interpolator
                itp = linear_interp_gpu(intt_in, bist_in)

                # loop over the length of the bisector
                for k in idz:sdz:CUDA.size(bisall_in, 1)
                    new_intt = (1.0 - dtrim) + (k-1) * step
                    if (1.0 - dtrim) >= CUDA.first(intt_in)
                        @inbounds bist_out[k] = itp(new_intt)
                    end
                    @inbounds intt_out[k] = new_intt
                    @inbounds widt_out[k] = widt_in[k]
                end
            end
        end
    end
    return nothing
end

function time_average_bis!(lenall, bisall_mean, intall_mean, widall_mean, bisall_in, intall_in, widall_in)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y
    # idz = threadIdx().z + blockDim().z * (blockIdx().z-1)
    # sdz = blockDim().z * gridDim().z

    # loop over disk positions for bisectors
    for i in idx:sdx:CUDA.length(lenall)

        # loop over intensity in bisector
        for k in idy:sdy:CUDA.size(bisall_in, 1)

            # holder for mean
            bis_sum = CUDA.zero(eltype(bisall_mean))
            int_sum = CUDA.zero(eltype(intall_mean))
            wid_sum = CUDA.zero(eltype(widall_mean))

            # loop over epochs of bisectors
            for j in 1:lenall[i]
                bis_sum += bisall_in[k, j, i]
                int_sum += intall_in[k, j, i]
                wid_sum += widall_in[k, j, i]
            end

            # take the mean and allocate 
            @inbounds bisall_mean[k, i] = bis_sum / lenall[i]
            @inbounds intall_mean[k, i] = int_sum / lenall[i]
            @inbounds widall_mean[k, i] = wid_sum / lenall[i]
        end        
    end
    return nothing 
end
