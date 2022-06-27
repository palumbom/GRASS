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

function sort_data_for_gpu(soldata::SolarData{T}) where T<:AbstractFloat
    # allocate memory for arrays to pass to gpu
    len = collect(values(soldata.len))
    npositions = length(len)
    wav = zeros(100, maximum(len), npositions)
    bis = zeros(100, maximum(len), npositions)
    wid = zeros(100, maximum(len), npositions)
    dep = zeros(100, maximum(len), npositions)
    for (ind,(key,val)) in enumerate(soldata.len)
        wav[:, 1:val, ind] .= soldata.wav[key]
        bis[:, 1:val, ind] .= soldata.bis[key]
        wid[:, 1:val, ind] .= soldata.wid[key]
        dep[:, 1:val, ind] .= soldata.dep[key]
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
    wav .= view(wav, :, :, inds_mu)
    bis .= view(bis, :, :, inds_mu)
    wid .= view(wid, :, :, inds_mu)
    dep .= view(dep, :, :, inds_mu)

    # get indices to sort by axis within mu sort
    for (idx, val) in enumerate(unique(disc_mu))
        inds1 = (disc_mu .== val)
        inds2 = sortperm(disc_ax[inds1])
        disc_mu[inds1] .= disc_mu[inds1][inds2]
        disc_ax[inds1] .= disc_ax[inds1][inds2]
    end
    return disc_mu, disc_ax, len, wav, bis, wid, dep
end

