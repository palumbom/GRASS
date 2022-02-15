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

    # calculate the data index
    if mu_ind == CUDA.length(disc_mu)
        return CUDA.length(disc_mu)
    else
        mu_orig = disc_mu[mu_ind]
        while (disc_ax[mu_ind] != ax_val)
            if mu > disc_mu[mu_ind]
                mu_ind -= 1
                if disc_mu[mu_ind] != mu_orig
                    return mu_ind + 1
                end
            else
                mu_ind += 1
                if disc_mu[mu_ind] != mu_orig
                    return mu_ind - 1
                end
            end
        end
        return mu_ind
    end
    return nothing
end

function sort_data_for_gpu(soldata::SolarData{T}) where T<:AbstractFloat
    # allocate memory for arrays to pass to gpu
    len = collect(values(soldata.len))
    wav = zeros(100, maximum(len), length(soldata.len))
    bis = zeros(100, maximum(len), length(soldata.len))
    wid = zeros(100, maximum(len), length(soldata.len))
    dep = zeros(100, maximum(len), length(soldata.len))
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
    wav .= wav[:, :, inds_mu]
    bis .= bis[:, :, inds_mu]
    wid .= wid[:, :, inds_mu]
    dep .= dep[:, :, inds_mu]

    # get indices to sort by axis within mu sort
    for (idx, val) in enumerate(unique(disc_mu))
        inds1 = (disc_mu .== val)
        inds2 = sortperm(disc_ax[inds1])
        disc_mu[inds1] .= disc_mu[inds1][inds2]
        disc_ax[inds1] .= disc_ax[inds1][inds2]
    end
    return disc_mu, disc_ax, len, wav, bis, wid, dep
end

