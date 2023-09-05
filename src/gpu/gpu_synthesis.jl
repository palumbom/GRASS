function fill_workspaces!(line, variability, extra_z, tloop, dat_idx, z_rot,
                          z_cbs, bisall, intall, widall, allwavs, allints)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(dat_idx)
        # move to next iter if off disk
        d_idx = dat_idx[i]
        if CUDA.iszero(d_idx)
            continue
        end

        # alias time index
        t = tloop[i]

        # calculate shifted line center
        λΔD = line * (1.0 + z_rot[i]) * (1.0 + z_cbs[i] * variability) * (1.0 + extra_z)

        # get length of input data arrays to loop over
        lent = 100
        for k in idy:sdy:lent
            # get forward and reverse indices
            idx1 = k
            idx2 = lent - (k - 1)

            # slice out the correct views of the input data for position
            @inbounds bis1 = bisall[idx1, t, d_idx]
            @inbounds wid1 = widall[idx1, t, d_idx]
            @inbounds int1 = intall[idx1, t, d_idx]

            @inbounds bis2 = bisall[idx2, t, d_idx]
            @inbounds wid2 = widall[idx2, t, d_idx]
            @inbounds int2 = intall[idx2, t, d_idx]

            # right side of line, indexing from middle left to right
            @inbounds allwavs[i, k+lent] = (λΔD + (0.5 * wid1 + bis1))
            @inbounds allints[i, k+lent] = int1

            # left sight of line, indexing from middle right to left
            @inbounds allwavs[i, k] = (λΔD - (0.5 * wid2 - bis2))
            @inbounds allints[i, k] = int2
        end
    end
    return nothing
end


function line_profile_gpu!(star_map, μs, lambdas, allwavs, allints)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(μs)
        # move to next iter if off disk
        if μs[i] <= 0.0
            continue
        end

        # take view of arrays to pass to interpolater
        allwavs_i = CUDA.view(allwavs, i, :)
        allints_i = CUDA.view(allints, i, :)

        # set up interpolator
        itp = linear_interp_gpu(allwavs_i, allints_i)

        # loop over wavelengths
        for k in idy:sdy:CUDA.length(lambdas)
            if ((lambdas[k] < CUDA.first(allwavs_i)) || (lambdas[k] > CUDA.last(allwavs_i)))
                continue
            else
                @inbounds star_map[i,k] *= itp(lambdas[k])
            end
        end
    end
    return nothing
end
