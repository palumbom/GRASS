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
    # get the mean of the first dataset
    group1 = view(arr, :, 1:ntimes[1])
    meangroup1 = dropdims(mean(group1, dims=2), dims=2)

    # loop over the nth datasets
    for i in 2:length(ntimes)
        # get the mena
        groupn = view(arr, :, sum(ntimes[1:i-1])+1:sum(ntimes[1:i]))
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
    bis .-= mean(bis, dims=2)
    return nothing
end
