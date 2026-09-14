# the _out buffers persist across lines, so every branch must write every element
function trim_bisector_gpu!(depth, variability, depcontrast, lenall, bisall_out,
                            intall_out, widall_out, bisall_in, intall_in, widall_in)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y
    idz = threadIdx().z + blockDim().z * (blockIdx().z-1)
    sdz = blockDim().z * gridDim().z

    # scalar arithmetic in the array precision: a Float64 depth or literal would
    # promote the kernel to double precision whatever precision was requested
    T = eltype(intall_in)

    # loop over disk positions for bisectors
    for i in idx:sdx:CUDA.length(lenall)
        # get depth to trim to
        dtrim = T(depth) * depcontrast[i]

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

            # chop resamples to maximum(intt), scale to 1.0; intt is ascending
            int_top = one(T)
            if (one(T) - dtrim) >= CUDA.first(intt_in)
                int_top = CUDA.last(intt_in)
            end
            step = (int_top - (one(T) - dtrim))/(CUDA.length(intt_in) - 1)

            if variability
                # set up interpolator
                itp = linear_interp_gpu(intt_in, bist_in)

                # loop over the length of the bisector
                for k in idz:sdz:CUDA.size(bisall_in, 1)
                    new_intt = (one(T) - dtrim) + (k-1) * step
                    if (one(T) - dtrim) >= CUDA.first(intt_in)
                        @inbounds bist_out[k] = itp(new_intt)
                    else
                        # scaling leaves the bisector untrimmed
                        @inbounds bist_out[k] = bist_in[k]
                    end
                    @inbounds intt_out[k] = new_intt

                    # variable widths pass through unchanged
                    @inbounds widall_out[k,j,i] = widall_in[k,j,i]
                end
            else
                for k in idz:sdz:CUDA.size(bisall_in, 1)
                    @inbounds intt_out[k] = (one(T) - dtrim) + (k-1) * step
                    @inbounds bist_out[k] = zero(T)
                    @inbounds widall_out[k,j,i] = widall_in[k,1,i]
                end
            end
        end
    end
    return nothing
end

# Collapse the time axis of the bisector input so every epoch of a tile carries that tile's
# time-averaged profile. The line keeps its asymmetry and its centre-to-limb variation but
# loses all granulation jitter, which isolates the two effects that `variability = true`
# otherwise turns on together.
#
# Operates on the raw input arrays before trimming, so everything downstream is unchanged: the
# tloop index still varies per cell, it just selects among identical epochs.
function broadcast_mean_bis!(lenall, bisall, intall, widall, bisall_mean, intall_mean, widall_mean)
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    for i in idx:sdx:CUDA.length(lenall)
        for k in idy:sdy:CUDA.size(bisall, 1)
            @inbounds b = bisall_mean[k, i]
            @inbounds v = intall_mean[k, i]
            @inbounds w = widall_mean[k, i]
            for j in 1:lenall[i]
                @inbounds bisall[k, j, i] = b
                @inbounds intall[k, j, i] = v
                @inbounds widall[k, j, i] = w
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
