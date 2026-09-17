# With interp or pool set, each cell's line shape is its nearest tile's current epoch plus the
# difference between a target time mean and the nearest tile's own time mean (`*_mean_own`).
# The target is the `*_mean_src` arrays interpolated over limb angle (tiles
# dat_idx_lo/dat_idx_hi, weight dat_wt) when interp is set, else the source at the nearest
# tile; with pool_axes the source is the axis-pooled mean, else the tile's own. The means are
# of the TRIMMED arrays, so every tile is on the same fractional-depth index and the arrays
# can be blended index by index; the blended intensity grid then carries the blended depth.
# Written as `q + (target - q_mean_k)` so that with no interpolation and no pooling the
# correction is exactly zero.
function fill_workspaces_2D_eclipse!(line, variability, extra_z, tloop, dat_idx, z_rot,
    z_cbs, lenall, bisall, intall, widall, allwavs, allints, contrast,
    interp, pool, dat_idx_lo, dat_idx_hi, dat_wt,
    bisall_mean_src, intall_mean_src, widall_mean_src,
    bisall_mean_own, intall_mean_own, widall_mean_own)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    T = eltype(bisall)
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

        # alias time index, wrapped into the data length of this tile's current key
        len = lenall[k]
        t = mod1(tloop[m,n], len)
        @inbounds tloop[m,n] = t

        # bracketing tiles and weight for the interpolated time-mean shape
        klo = k
        khi = k
        w = zero(T)
        if interp
            @inbounds klo = dat_idx_lo[m,n]
            @inbounds khi = dat_idx_hi[m,n]
            @inbounds w = dat_wt[m,n]
        end

        # absolute wavelengths are computed in Float64 on purpose; storing them in a
        # Float32 allwavs is the precision-limiting step of the single-precision path
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

            if interp | pool
                @inbounds bis1 += ((one(T) - w) * bisall_mean_src[idx1, klo] + w * bisall_mean_src[idx1, khi]) - bisall_mean_own[idx1, k]
                @inbounds wid1 += ((one(T) - w) * widall_mean_src[idx1, klo] + w * widall_mean_src[idx1, khi]) - widall_mean_own[idx1, k]
                @inbounds int1 += ((one(T) - w) * intall_mean_src[idx1, klo] + w * intall_mean_src[idx1, khi]) - intall_mean_own[idx1, k]
                @inbounds bis2 += ((one(T) - w) * bisall_mean_src[idx2, klo] + w * bisall_mean_src[idx2, khi]) - bisall_mean_own[idx2, k]
                @inbounds wid2 += ((one(T) - w) * widall_mean_src[idx2, klo] + w * widall_mean_src[idx2, khi]) - widall_mean_own[idx2, k]
                @inbounds int2 += ((one(T) - w) * intall_mean_src[idx2, klo] + w * intall_mean_src[idx2, khi]) - intall_mean_own[idx2, k]
            end

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