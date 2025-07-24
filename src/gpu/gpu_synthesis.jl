function fill_workspaces!(line, variability, extra_z, tloop, dat_idx, z_rot,
                          z_cbs, lenall, bisall, intall, widall, allwavs, allints)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(dat_idx)
        # move to next iter if off disk
        k = dat_idx[i]
        if CUDA.iszero(k)
            continue
        end

        # alias time index
        len = lenall[dat_idx[i]]
        if tloop[i] > len
            @inbounds tloop[i] -= len
        end
        t = tloop[i]

        # calculate shifted line center
        λΔD = line * (1.0 + z_rot[i]) * (1.0 + z_cbs[i] * variability) * (1.0 + extra_z)

        # get length of input data arrays to loop over
        lent = 100
        for j in idy:sdy:lent
            # get forward and reverse indices
            idx1 = j
            idx2 = lent - (j - 1)

            # slice out the correct views of the input data for position
            @inbounds bis1 = bisall[idx1, t, k]
            @inbounds wid1 = widall[idx1, t, k]
            @inbounds int1 = intall[idx1, t, k]

            @inbounds bis2 = bisall[idx2, t, k]
            @inbounds wid2 = widall[idx2, t, k]
            @inbounds int2 = intall[idx2, t, k]

            # right side of line, indexing from middle left to right
            @inbounds allwavs[i, j+lent] = (λΔD + (0.5 * wid1 + bis1))
            @inbounds allints[i, j+lent] = int1

            # left sight of line, indexing from middle right to left
            @inbounds allwavs[i, j] = (λΔD - (0.5 * wid2 - bis2))
            @inbounds allints[i, j] = int2
        end
    end
    return nothing
end


function line_profile_gpu!(prof, μs, wts, λs, allwavs, allints)
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
        for j in idy:sdy:CUDA.length(λs)
            if ((λs[j] < CUDA.first(allwavs_i)) || (λs[j] > CUDA.last(allwavs_i)))
                @inbounds CUDA.@atomic prof[j] += wts[i]
            else
                @inbounds CUDA.@atomic prof[j] += itp(λs[j]) * wts[i]
            end
        end
    end
    return nothing
end

function line_profile_resolved_gpu!(prof, μ_bins, μs, wts, λs, allwavs, allints)
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

        # get mu_idx
        μ_idx = searchsortednearest_gpu(μ_bins, μs[i])

        # take view of arrays to pass to interpolater
        allwavs_i = CUDA.view(allwavs, i, :)
        allints_i = CUDA.view(allints, i, :)

        # set up interpolator
        itp = linear_interp_gpu(allwavs_i, allints_i)

        # loop over wavelengths
        for j in idy:sdy:CUDA.length(λs)
            if ((λs[j] < CUDA.first(allwavs_i)) || (λs[j] > CUDA.last(allwavs_i)))
                @inbounds CUDA.@atomic prof[j, μ_idx] += wts[i]
            else
                @inbounds CUDA.@atomic prof[j, μ_idx] += itp(λs[j]) * wts[i]
            end
        end
    end
    return nothing
end

function apply_line!(t, prof, flux, sum_wts)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(prof)
        @inbounds flux[i,t] *= prof[i] / sum_wts
        @inbounds prof[i] = 0.0
    end
    return nothing
end

function apply_line_resolved!(t, prof, flux, sum_wts)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.size(prof,1)
        for j in idy:sdy:CUDA.size(prof,2)
            @inbounds flux[i,j,t] *= prof[i,j] / sum_wts
            @inbounds prof[i,j] = 0.0
        end
    end
    return nothing
end
