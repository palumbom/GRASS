function sort_data_for_gpu(soldata::SolarData{T}) where T<:AbstractFloat
    # allocate memory for arrays to pass to gpu
    len = collect(values(soldata.len))
    cbs = collect(values(soldata.cbs))
    npositions = length(len)
    bis = zeros(100, maximum(len), npositions)
    int = zeros(100, maximum(len), npositions)
    wid = zeros(100, maximum(len), npositions)
    for (ind,(key,val)) in enumerate(soldata.len)
        bis[:, 1:val, ind] .= soldata.bis[key]
        int[:, 1:val, ind] .= soldata.int[key]
        wid[:, 1:val, ind] .= soldata.wid[key]
    end

    # get the value of mu and ax codes
    disc_ax = parse_ax_string.(getindex.(keys(soldata.len),1))
    disc_mu = parse_mu_string.(getindex.(keys(soldata.len),2))

    # get indices to sort by mus
    inds_mu = sortperm(disc_mu)
    disc_mu .= disc_mu[inds_mu]
    disc_ax .= disc_ax[inds_mu]

    # get the arrays in mu sorted order
    len .= len[inds_mu]
    cbs .= cbs[inds_mu]
    bis .= view(bis, :, :, inds_mu)
    int .= view(int, :, :, inds_mu)
    wid .= view(wid, :, :, inds_mu)

    # get indices to sort by axis within mu sort
    for mu_val in unique(disc_mu)
        inds1 = (disc_mu .== mu_val)
        inds2 = sortperm(disc_ax[inds1])
        disc_mu[inds1] .= disc_mu[inds1][inds2]
        disc_ax[inds1] .= disc_ax[inds1][inds2]

        len[inds1] .= len[inds1][inds2]
        cbs[inds1] .= cbs[inds1][inds2]
        bis[:, :, inds1] .= bis[:, :, inds1][:, :, inds2]
        int[:, :, inds1] .= int[:, :, inds1][:, :, inds2]
        wid[:, :, inds1] .= wid[:, :, inds1][:, :, inds2]
    end
    return disc_mu, disc_ax, len, cbs, bis, int, wid
end

function find_nearest_ax_gpu(x::T, y::T) where T<:AbstractFloat
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

function find_data_index_gpu(x, y, disc_mu, disc_ax)
    # find the nearest mu ind and ax code
    mu = calc_mu(x,y)
    mu_ind = searchsortednearest_gpu(disc_mu, mu)
    ax_val = find_nearest_ax_gpu(x, y)

    # find the first index of disc_mu with that discrete mu val
    i = 1
    while disc_mu[i] != disc_mu[mu_ind]
        i += 1
    end
    mu_ind = i

    # calculate the data index
    if mu_ind == CUDA.length(disc_mu)
        # return immediately if nearest mu is disk center
        return CUDA.length(disc_mu)
    else
        # otherwise we need the right axis value
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
    end
    return nothing
end

function iterate_tloop_gpu!(tloop, data_inds, lenall, grid)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            # find position on disk and move to next iter if off disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
            if r2 > 1.0
                continue
            end

            # check that tloop didn't overshoot the data and iterate
            ntimes = lenall[data_inds[i,j]]
            if tloop[i,j] < ntimes
                @inbounds tloop[i,j] += 1
            else
                @inbounds tloop[i,j] = 1
            end
        end
    end
    return nothing
end

function initialize_arrays_for_gpu(data_inds, tloop, norm_terms, z_rot, z_cbs,
                                   grid, disc_mu, disc_ax, lenall, cbsall,
                                   u1, u2, polex, poley, polez)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # make some aliases
    len = CUDA.length(grid)
    rstar = one(eltype(grid))

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            # find position on disk and move to next iter if off disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
            if r2 > 1.0
                continue
            end

            # find the correct data index and initialize tloop value
            idx = find_data_index_gpu(x, y, disc_mu, disc_ax)
            @inbounds data_inds[i,j] = idx
            @inbounds tloop[i,j] = CUDA.floor(Int32, rand() * lenall[idx]) + 1

            # calculate the normalization
            @inbounds norm_terms[i,j] = calc_norm_term(x, y, len, u1, u2)

            # calculate the rotational and convective doppler shift
            @inbounds z_rot[i,j] = patch_velocity_los_gpu(x, y, rstar, polex, poley, polez)
            @inbounds z_cbs[i,j] = cbsall[idx]
        end
    end
    return nothing
end
