function trim_bisector_gpu!(depth, bisall_out, intall_out, bisall_in, intall_in)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # loop over epochs of bisectors
    for i in idx:sdx:CUDA.size(bisall_in, 2)
        # get views of the correct time slice in input
        bist_in = CUDA.view(bisall_in, :, i)
        intt_in = CUDA.view(intall_in, :, i)

        # get views of the correct time slice in output
        bist_out = CUDA.view(bisall_out, :, i)
        intt_out = CUDA.view(intall_out, :, i)

        # loop over the length of the bisector
        step = depth/(CUDA.length(intt_in) - 1)
        for j in idy:sdy:CUDA.size(bisall_in, 1)
            # set the new depth value
            new_intt = (1.0 - depth) + (j-1) * step

            # if depth is greater than input depth, stretch the bisector
            if (1.0 - depth) < CUDA.first(intt_in)
                @inbounds intt_out[j] = new_intt
            else
                # interpolate to get the new wavelength value
                if new_intt <= CUDA.first(intt_in)
                    @inbounds bist_out[j] = CUDA.first(bist_in)
                elseif new_intt >= CUDA.last(intt_in)
                    @inbounds bist_out[j] = CUDA.last(bist_in)
                else
                    k = CUDA.searchsortedfirst(intt_in, new_intt) - 1
                    k0 = CUDA.clamp(k, CUDA.firstindex(bist_in), CUDA.lastindex(bist_in))
                    k1 = CUDA.clamp(k+1, CUDA.firstindex(bist_in), CUDA.lastindex(bist_in))
                    @inbounds bist_out[j] = (bist_in[k0] * (intt_in[k1] - new_intt) + bist_in[k1] * (new_intt - intt_in[k0])) / (intt_in[k1] - intt_in[k0])
                end

                # assign bisector fluxes from dept
                @inbounds intt_out[j] = new_intt
            end
        end
    end
    return nothing
end
