function find_nearest_ax_gpu(x, y)
    if (CUDA.iszero(x) & CUDA.iszero(y))
        return 0 # center
    elseif y >= CUDA.abs(x)
        return 1 # north
    elseif y <= -CUDA.abs(x)
        return 2 # south
    elseif x <= -CUDA.abs(y)
        return 3 # east
    elseif x >= CUDA.abs(y)
        return 4 # west
    else
        return 0
    end
end

function find_data_index_gpu(μ, ax_val, disc_mu, disc_ax)
    # find the nearest mu ind and ax code
    mu_ind = searchsortednearest_gpu(disc_mu, μ)

    # return immediately if nearest mu is disk center
    if mu_ind == CUDA.length(disc_mu)
        return CUDA.length(disc_mu)
    end

    # find the first index of disc_mu with that discrete mu val
    i = 1
    while disc_mu[i] != disc_mu[mu_ind]
        i += 1
    end
    mu_ind = i

    # get the right axis value
    mu_ind_orig = mu_ind
    mu_val_orig = disc_mu[mu_ind_orig]
    while ((disc_ax[mu_ind] != ax_val) & (disc_mu[mu_ind] == mu_val_orig))
        mu_ind += 1
    end

    # check that we haven't overflowed into the next batch of mus
    if disc_mu[mu_ind] == mu_val_orig
        return mu_ind
    else
        return mu_ind_orig
    end
    return nothing
end

function iterate_tloop_gpu!(tloop, dat_idx, lenall)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x

    # parallelized loop over grid
    for i in idx:sdx:CUDA.size(dat_idx,1)
        if CUDA.iszero(dat_idx[i])
            continue
        end

        # check that tloop didn't overshoot the data and iterate
        ntimes = lenall[dat_idx[i]]
        if tloop[i] < ntimes
            @inbounds tloop[i] += 1
        else
            @inbounds tloop[i] = 1
        end
    end
    return nothing
end

function check_tloop_gpu!(tloop, dat_idx, lenall)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(dat_idx,)
        if CUDA.iszero(dat_idx[i])
            continue
        end

        # check that tloop didn't overshoot the data and iterate
        ntimes = lenall[dat_idx[i]]
        if tloop[i] > ntimes
            @inbounds tloop[i] = 1
        end
    end
    return nothing
end

function generate_tloop_gpu!(tloop::AA{Int32,1}, gpu_allocs::GPUAllocs{T}, soldata::GPUSolarData{T}) where T<:AF
    dat_idx = gpu_allocs.dat_idx
    lenall = soldata.len

    threads1 = 256
    blocks1 = cld(CUDA.length(dat_idx), prod(threads1))

    @cusync @captured @cuda threads=threads1 blocks=blocks1 generate_tloop_gpu!(tloop, dat_idx, lenall)
    return nothing
end

function generate_tloop_gpu!(tloop, dat_idx, lenall)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(dat_idx)
        if CUDA.iszero(dat_idx[i])
            continue
        end
        idx = dat_idx[i]
        @inbounds tloop[i] = CUDA.floor(Int32, rand() * lenall[idx]) + 1
    end
    return nothing
end
