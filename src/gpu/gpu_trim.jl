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

        # set up interpolator
        itp = linear_interp_gpu(intt_in, bist_in)

        # get step size for loop
        step = depth/(CUDA.length(intt_in) - 1)

        # loop over the length of the bisector
        for j in idy:sdy:CUDA.size(bisall_in, 1)
            new_intt = (1.0 - depth) + (j-1) * step
            if (1.0 - depth) >= CUDA.first(intt_in)
                @inbounds bist_out[j] = itp(new_intt)
            end
            @inbounds intt_out[j] = new_intt
        end
    end
    return nothing
end
