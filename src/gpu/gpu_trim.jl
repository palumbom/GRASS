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

            # Resampled intensity grid endpoint is branch-dependent, mirroring the CPU:
            # trim_bisector_chop! runs to maximum(intt), trim_bisector_scale! to 1.0.
            # They coincide only for templates whose intensities reach the continuum;
            # FeI_5382 tops out at 0.985, and using 1.0 there costs 44 m/s against the CPU.
            # intt is ascending (linear_interp_gpu below already requires it), so the
            # maximum is the last element -- do not replace this with a scan, every
            # z-thread would repeat it.
            int_top = 1.0
            if (1.0 - dtrim) >= CUDA.first(intt_in)
                int_top = CUDA.last(intt_in)
            end
            step = (int_top - (1.0 - dtrim))/(CUDA.length(intt_in) - 1)

            if variability
                # set up interpolator
                itp = linear_interp_gpu(intt_in, bist_in)

                # loop over the length of the bisector
                for k in idz:sdz:CUDA.size(bisall_in, 1)
                    new_intt = (1.0 - dtrim) + (k-1) * step
                    if (1.0 - dtrim) >= CUDA.first(intt_in)
                        @inbounds bist_out[k] = itp(new_intt)
                    else
                        # scaling, not chopping: trim_bisector_scale! rewrites only the
                        # intensities, leaving the bisector at its untrimmed input value.
                        # bisall_out persists across lines, so this must be written
                        # explicitly or line l inherits line l-1's chopped bisector.
                        @inbounds bist_out[k] = bist_in[k]
                    end
                    @inbounds intt_out[k] = new_intt

                    # widall_out has the same cross-line lifetime as bisall_out, and the
                    # non-variability branch below overwrites every epoch with epoch 1.
                    # A variable line must restore its own epoch's widths or it inherits
                    # the fixed ones -- across lines and across time steps.
                    @inbounds widall_out[k,j,i] = widall_in[k,j,i]
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
