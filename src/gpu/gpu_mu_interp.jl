# Tile lookup for interpolating the template data over limb angle. For every on-disk cell,
# find the two tiles whose mu bracket the cell's mu, each on the cell's own axis when that axis
# exists at that mu level and otherwise on the level's fallback tile (the same fallback
# find_data_index_gpu uses), and the weight of the upper tile, linear in theta = acos(mu).
# Cells beyond either end of the mu grid are clamped to the end tile with weight zero: there
# is no extrapolation below the lowest tile. z_cbs is overwritten with the interpolated
# convective blueshift.
#
# The nearest tile in dat_idx is left alone. The granulation draw (tloop) and the per-epoch
# departure from each tile's time mean still come from that tile; only the time-mean line
# shape is interpolated (fill_workspaces_2D_eclipse!).
# Axis-pooled time means: every tile at a mu level receives the unweighted mean over that
# level's tiles of the `own` means, so the shape a cell is corrected toward is axisymmetric
# while its granulation departure still comes from its own tile (pool_axes).
function pool_axis_means_gpu!(bis_pool, int_pool, wid_pool, bis_own, int_own, wid_own, disc_mu)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    T = eltype(bis_pool)

    for k in idx:sdx:CUDA.size(bis_own, 2)
        for j in idy:sdy:CUDA.size(bis_own, 1)
            bis_sum = zero(T)
            int_sum = zero(T)
            wid_sum = zero(T)
            count = 0
            for i in 1:CUDA.size(bis_own, 2)
                if disc_mu[i] == disc_mu[k]
                    @inbounds bis_sum += bis_own[j, i]
                    @inbounds int_sum += int_own[j, i]
                    @inbounds wid_sum += wid_own[j, i]
                    count += 1
                end
            end
            @inbounds bis_pool[j, k] = bis_sum / count
            @inbounds int_pool[j, k] = int_sum / count
            @inbounds wid_pool[j, k] = wid_sum / count
        end
    end
    return nothing
end

function get_interp_keys_gpu!(dat_idx_lo, dat_idx_hi, dat_wt, z_cbs, μs, ax_codes, cbsall, disc_mu, disc_ax)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x

    T = eltype(dat_wt)
    n = CUDA.length(disc_mu)

    for i in idx:sdx:CUDA.length(μs)
        μ = μs[i]

        # off-disk cells get no tiles, matching dat_idx
        if μ <= zero(T)
            @inbounds dat_idx_lo[i] = 0
            @inbounds dat_idx_hi[i] = 0
            @inbounds dat_wt[i] = zero(T)
            continue
        end

        if μ >= disc_mu[n]
            # at or above disk centre, which is the single last tile
            lo = n
            hi = n
            w = zero(T)
        elseif μ <= disc_mu[1]
            # below the lowest tile: clamp, do not extrapolate
            lo = find_data_index_at_mu_gpu(1, ax_codes[i], disc_mu, disc_ax)
            hi = lo
            w = zero(T)
        else
            # disc_mu is sorted ascending with repeats, so ih is the first tile of the level
            # at or above μ and ih - 1 the last tile of the level strictly below it
            ih = CUDA.searchsortedfirst(disc_mu, μ)
            il = ih - 1
            θ_lo = CUDA.acos(disc_mu[il])
            θ_hi = CUDA.acos(disc_mu[ih])
            w = (θ_lo - CUDA.acos(μ)) / (θ_lo - θ_hi)
            lo = find_data_index_at_mu_gpu(il, ax_codes[i], disc_mu, disc_ax)
            hi = find_data_index_at_mu_gpu(ih, ax_codes[i], disc_mu, disc_ax)
        end

        @inbounds dat_idx_lo[i] = lo
        @inbounds dat_idx_hi[i] = hi
        @inbounds dat_wt[i] = w
        @inbounds z_cbs[i] = (one(T) - w) * cbsall[lo] + w * cbsall[hi]
    end
    return nothing
end
