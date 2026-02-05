function fill_workspaces_2D_eclipse!(line, variability, extra_z, tloop, dat_idx, z_rot,
    z_cbs, lenall, bisall, intall, widall, allwavs, allints, contrast)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    Nθ_max = CUDA.size(dat_idx, 2)

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(dat_idx)
        row = (i - 1) ÷ Nθ_max
        col = (i - 1) % Nθ_max
        m = row + 1
        n = col + 1

        # move to next iter if off disk
        k = dat_idx[m,n]
        if CUDA.iszero(k)
            continue
        end

        # alias time index
        len = lenall[k]
        if tloop[m,n] > len
            @inbounds tloop[m,n] -= len
        end
        t = tloop[m,n]

        # calculate shifted line center
        λΔD = line * (1.0 + z_rot[m,n]) * (1.0 + z_cbs[m,n] * variability * contrast[m,n]) * (1.0 + extra_z * variability * contrast[m,n])

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
            @inbounds allwavs[m, n, j+lent] = (λΔD + (0.5 * wid1 + bis1))
            @inbounds allints[m, n, j+lent] = int1

            # left sight of line, indexing from middle right to left
            @inbounds allwavs[m, n, j] = (λΔD - (0.5 * wid2 - bis2))
            @inbounds allints[m, n, j] = int2
        end
    end
    return nothing
end