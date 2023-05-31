function parse_mu_string(s::String)
    s = s[3:end]
    return tryparse(Float64, s[1] * "." * s[2:end])
end

function mu_to_symb(mu::T) where T<:AF
    return Symbol("mu" * replace(string(mu), "."=>""))
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

function adjust_data_mean(arr::AA{T,2}, ntimes::Vector{Int64};
                          lo_ind::Int=15, hi_ind::Int=75) where T<:Real
    # get indices for number of contiguous obs
    arr_idx = vcat([0], cumsum(ntimes))

    # find the mean of the longest data set, use bottom n% of bisector only
    idx = argmax(ntimes)
    arr_sub = view(arr, lo_ind:hi_ind, arr_idx[idx]+1:arr_idx[idx+1])
    mean_ref = mean(arr_sub)

    # loop over the datasets and adjust mean
    for i in 1:length(ntimes)
        # get the group of measurements
        group = view(arr, :, arr_idx[i]+1:arr_idx[i+1])

        # get the mean of bottom n% of bisctor
        mean_group = mean(view(group, lo_ind:hi_ind, :))

        # find the distance between the means and correct by it
        mean_dist = mean_ref - mean_group
        group .+= mean_dist
    end
    return nothing
end

function identify_bad_cols(bisall::AA{T,2}, intall::AA{T,2}, widall::AA{T,2};
                           lo_ind::Int=15, hi_ind::Int=75) where T<:AF
    @assert size(bisall) == size(intall) == size(widall)

    # how many sigma away to consider outlier
    nsigma = 2.0

    # allocate boolean array (column will be stripped if ool[i] == true)
    badcols = zeros(Bool, size(bisall,2))

    # make sure the max/min width is reasonable
    max_wid_view = view(widall, size(widall, 1), :)
    max_wid_avg = mean(max_wid_view)
    max_wid_std = std(max_wid_view)

    min_wid_view = view(widall, 1, :)
    min_wid_avg = mean(min_wid_view)
    min_wid_std = std(min_wid_view)

    # remove significant max width outliers
    idx1 = abs.(max_wid_avg .- max_wid_view) .> (nsigma .* max_wid_std)
    idx2 = abs.(min_wid_avg .- min_wid_view) .> (nsigma .* min_wid_std)
    badcols[idx1] .= true

    # get views
    bis_view = view(bisall, lo_ind:hi_ind, .!badcols)
    wid_view = view(widall, lo_ind:hi_ind, .!badcols)

    # take averages
    bis_avg = dropdims(mean(bis_view, dims=2), dims=2)
    wid_avg = dropdims(mean(wid_view, dims=2), dims=2)
    bis_std = dropdims(std(bis_view, dims=2), dims=2)
    wid_std = dropdims(std(wid_view, dims=2), dims=2)

    # loop through checking for bad columns
    for i in 1:size(bisall,2)
        # if column is already marked bad, move on
        if badcols[i]
            continue
        end

        # get views of data
        bist = view(bisall, lo_ind:hi_ind, i)
        intt = view(intall, lo_ind:hi_ind, i)
        widt = view(widall, lo_ind:hi_ind, i)

        # check for monotinicity in measurements
        if !ismonotonic(widt) | !ismonotonic(intt)
            badcols[i] = true
        end

        # check for skipped epochs in preprocessing
        if all(iszero.(intt))
            badcols[i] = true
        end

        # check for NaNs
        if any(isnan.(intt))
            badcols[i] = true
        end

        # remove measurements that are significant outliers
        bis_cond = any(abs.(bis_avg .- bist) .> (nsigma .* bis_std))
        wid_cond = any(abs.(wid_avg .- widt) .> (nsigma .* wid_std))
        if bis_cond | wid_cond
            badcols[i] = true
        end
    end
    return badcols
end

function relative_bisector_wavelengths(bis::AA{T,2}) where T<:AF
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
    idx2 = searchsortedfirst(intt, top) - 2
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
