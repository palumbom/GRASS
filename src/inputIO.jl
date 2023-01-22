function parse_mu_string(s::String)
    s = s[3:end]
    return tryparse(Float64, s[1] * "." * s[2:end])
end

function parse_mu_string(s::Symbol)
    return parse_mu_string(string(s))
end

function parse_ax_string(s::String)
    if s == "c"; return 0; end;
    if s == "n"; return 1; end;
    if s == "s"; return 2; end;
    if s == "e"; return 3; end;
    if s == "w"; return 4; end;
end

function parse_ax_string(s::Symbol)
    return parse_ax_string(string(s))
end

function adjust_data_mean(arr::AA{T,2}, ntimes::Vector{Int64}) where T<:Real
    # get the mean of all the data
    mean_all = mean(arr)
    # mean_all = dropdims(mean(view(arr, :, 1:ntimes[1]), dims=2), dims=2)

    # loop over the nth datasets
    arr_idx = vcat([0], cumsum(ntimes))
    for i in 1:length(ntimes)
        # get the group of measurements
        group = view(arr, :, arr_idx[i]+1:arr_idx[i+1])

        # get the mean of the group
        mean_group = mean(group)
        # mean_group = dropdims(mean(group, dims=2), dims=2)

        # find the distance between the means and correct by it
        mean_dist = mean_all - mean_group
        group .-= mean_dist
    end
    return nothing
end

function adjust_data_mean_old(arr::AA{T,2}, ntimes::Vector{Int64}) where T<:Real
    # get the mean of the first dataset
    group1 = view(arr, :, 1:ntimes[1])
    meangroup1 = dropdims(mean(group1, dims=2), dims=2)

    # loop over the nth datasets
    arr_idx = vcat([0], cumsum(ntimes))
    for i in 2:length(ntimes)
        # get the mena
        groupn = view(arr, :, arr_idx[i]+1:arr_idx[i+1])
        meangroupn = dropdims(mean(groupn, dims=2), dims=2)

        # find the distance between the means and correct by it
        meandist = meangroupn - meangroup1
        groupn .-= meandist
    end
    return nothing
end

function identify_bad_cols(bisall::AA{T,2}, intall::AA{T,2}, widall::AA{T,2}) where T<:AF
    @assert size(bisall) == size(intall) == size(widall)

    # make boolean array (column will be stripped if badcol[i] == true)
    badcols = zeros(Bool, size(bisall,2))

    # find standarad deviation of data
    bis_std = dropdims(std(bisall, dims=2), dims=2)
    wid_std = dropdims(std(widall, dims=2), dims=2)

    # find mean and median of data
    bis_avg = dropdims(mean(bisall, dims=2), dims=2)
    wid_avg = dropdims(mean(widall, dims=2), dims=2)
    bis_med = dropdims(median(bisall, dims=2), dims=2)
    wid_med = dropdims(median(widall, dims=2), dims=2)

    # loop through checking for bad columns
    for i in 1:size(bisall,2)
        bist = view(bisall, :, i)
        intt = view(intall, :, i)
        widt = view(widall, :, i)

        # check for monotinicity in measurements
        if !ismonotonic(widt)
            badcols[i] = true
        end

        if !ismonotonic(intt)
            badcols[i] = true
        end

        # check for skipped epochs in preprocessing
        if all(iszero.(intt))
            badcols[i] = true
        end

        # remove data that is significant outlier
        nsigma = 3.0
        bis_cond = any(abs.(bis_avg .- bist) .> (nsigma .* bis_std))
        wid_cond = any(abs.(wid_avg .- widt) .> (nsigma .* wid_std))
        if bis_cond | wid_cond
            badcols[i] = true
        end
    end
    return badcols
end

function relative_bisector_wavelengths(bis::AA{T,2}) where T<:AF
    # TODO: find better way; how does this jive w/ convective blueshift?
    bis .-= mean(bis)
    return nothing
end

function extrapolate_input_data(bist::AA{T,1}, intt::AA{T,1}, widt::AA{T,1}, top::T) where T<:AF
    # fit the bottom bisector area and replace with model fit
    bot = minimum(intt)
    dep = 1.0 - bot
    idx1 = searchsortedfirst(intt, bot + 0.05 * dep)
    idx2 = searchsortedfirst(intt, bot + 0.15 * dep)
    bfit = pfit(view(intt, idx1:idx2), view(bist, idx1:idx2), 1)
    bist[1:idx1] .= bfit.(view(intt, 1:idx1))

    # fit the top bisector area and replace with model fit
    idx1 = searchsortedfirst(intt, top - 0.1 * dep )
    idx2 = searchsortedfirst(intt, top) - 1
    bfit = pfit(view(intt, idx1:idx2), view(bist, idx1:idx2), 1)
    bist[idx2:end] .= bfit.(view(intt, idx2:length(intt)))

    # extrapolate the width up to the continuum
    # TODO revisit this
    idx1 = length(widt) - 1
    wfit = pfit(view(intt, idx1:length(intt)), view(widt, idx1:length(intt)), 1)
    widt[idx1:end] .= wfit.([intt[idx1], 1.0])
    intt[end] = 1.0
    return nothing
end

function extrapolate_input_data(bis::AA{T,2}, int::AA{T,2}, wid::AA{T,2}, top::AA{T,1}) where T<:AF
    for i in 1:size(bis,2)
        # take a slice for one time snapshot
        bist = view(bis, :, i)
        intt = view(int, :, i)
        widt = view(wid, :, i)

        extrapolate_input_data(bist, intt, widt, top[i])
    end
    return nothing
end
