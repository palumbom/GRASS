# sort data for the GPU
function get_disc_mu_for_gpu(soldata::SolarData{T}; to_gpu::Bool=false) where T<:Float64
    mu_strings = string.(getindex.(keys(soldata.len),2))
    mu_strings .= [s[3:end][1] .* "." .* s[3:end][2:end] for s in mu_strings]
    disc_mu = parse.(Float64, mu_strings)
    if to_gpu
        return CuArray(disc_mu)
    end
    return disc_mu
end

function get_disc_ax_for_gpu(soldata::SolarData{T}; to_gpu::Bool=false) where T<:Float64
    disc_ax = string.(getindex.(keys(soldata.len),1))
    function parse_ax_strings(s::String)
        if s == "c"; return 0; end;
        if s == "n"; return 1; end;
        if s == "s"; return 2; end;
        if s == "e"; return 3; end;
        if s == "w"; return 4; end;
    end
    disc_ax = parse_ax_strings.(disc_ax)
    if to_gpu
        return CuArray(disc_ax)
    end
    return disc_ax
end

function find_nearest_ax_gpu(x::T, y::T) where T<:Float64
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

function find_data_index_gpu(mu_ind, ax_code)
    if mu_ind == 41
        return 41
    elseif !CUDA.iszero(CUDA.mod(mu_ind - 1, 4))
        return mu_ind + ax_code - 4
    else
        return mu_ind + ax_code - 1
    end
end

function sort_data_for_gpu(soldata::SolarData{T}) where T<:Float64
    # allocate memory for arrays to pass to gpu
    len = collect(values(soldata.len))
    wav = zeros(100, maximum(len), 41)
    bis = zeros(100, maximum(len), 41)
    wid = zeros(100, maximum(len), 41)
    dep = zeros(100, maximum(len), 41)
    for (ind,(key,val)) in enumerate(soldata.len)
        wav[:, 1:val, ind] .= soldata.wav[key]
        bis[:, 1:val, ind] .= soldata.bis[key]
        wid[:, 1:val, ind] .= soldata.wid[key]
        dep[:, 1:val, ind] .= soldata.dep[key]
    end

    # get the value of mu and ax codes
    disc_mu = get_disc_mu_for_gpu(soldata, to_gpu=false)
    disc_ax = get_disc_ax_for_gpu(soldata, to_gpu=false)

    # get keys to help debug
    the_keys = collect(keys(soldata.wav))

    # sort by mus
    inds_mu = sortperm(disc_mu)
    disc_mu .= disc_mu[inds_mu]
    disc_ax .= disc_ax[inds_mu]
    the_keys .= the_keys[inds_mu]

    # get the arrays in mu sorted order
    len .= len[inds_mu]
    wav .= wav[:, :, inds_mu]
    bis .= bis[:, :, inds_mu]
    wid .= wid[:, :, inds_mu]
    dep .= dep[:, :, inds_mu]

    # now sort by axis within mu sort
    for i in 0:length(unique(disc_mu)[1:end-1])-1
        # get inds to sort on
        inds_ax = sortperm(disc_ax[4i+1:4*(i+1)])

        # sort the disc arrays
        disc_mu[4i+1:4*(i+1)] .= disc_mu[4i+1:4*(i+1)][inds_ax]
        disc_ax[4i+1:4*(i+1)] .= disc_ax[4i+1:4*(i+1)][inds_ax]
        the_keys[4i+1:4*(i+1)] .= the_keys[4i+1:4*(i+1)][inds_ax]

        # sort the input data
        len[4i+1:4*(i+1)] .= len[4i+1:4*(i+1)][inds_ax]
        wav[:,:,4i+1:4*(i+1)] .= wav[:,:,4i+1:4*(i+1)][:, :, inds_ax]
        bis[:,:,4i+1:4*(i+1)] .= bis[:,:,4i+1:4*(i+1)][:, :, inds_ax]
        wid[:,:,4i+1:4*(i+1)] .= wid[:,:,4i+1:4*(i+1)][:, :, inds_ax]
        dep[:,:,4i+1:4*(i+1)] .= dep[:,:,4i+1:4*(i+1)][:, :, inds_ax]
    end
    return disc_mu, disc_ax, len, wav, bis, wid, dep
end

