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

function sort_mu_and_ax(soldata::SolarData{T}) where T<:AF
    # get the value of mu and ax codes
    disc_ax = parse_ax_string.(getindex.(keys(soldata.len),1))
    disc_mu = parse_mu_string.(getindex.(keys(soldata.len),2))

    # get indices to sort by mus
    inds_mu = sortperm(disc_mu)
    disc_mu .= disc_mu[inds_mu]
    disc_ax .= disc_ax[inds_mu]

    # get indices to sort by axis within mu sort
    for mu_val in unique(disc_mu)
        inds1 = (disc_mu .== mu_val)
        inds2 = sortperm(disc_ax[inds1])

        disc_mu[inds1] .= disc_mu[inds1][inds2]
        disc_ax[inds1] .= disc_ax[inds1][inds2]
    end
    return disc_mu, disc_ax
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

function identify_bad_cols(bisall::AA{T,2}, intall1::AA{T,2},
                           widall::AA{T,2}, intall2::AA{T,2};
                           lo_ind::Int=15, hi_ind::Int=75) where T<:AF
    @assert size(bisall) == size(intall1) == size(widall) == size(intall2)

    # how many sigma away to consider outlier
    nsigma = 4.0

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
        widt = view(widall, lo_ind:hi_ind, i)
        intt1 = view(intall1, lo_ind:hi_ind, i)
        intt2 = view(intall2, lo_ind:hi_ind, i)

        # check for monotinicity in measurements
        if !ismonotonic(widt) | !ismonotonic(intt1) | !ismonotonic(intt2)
            badcols[i] = true
        end

        # check for skipped epochs in preprocessing
        if all(iszero.(intt1)) | all(iszero.(intt2))
            badcols[i] = true
        end

        # check for NaNs
        if any(isnan.(intt1)) | any(isnan.(intt2))
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

function extrapolate_input_data(bist::AA{T,1}, intt1::AA{T,1},
                                widt::AA{T,1}, intt2::AA{T,1},
                                mu::T; weights=ones(length(bist))) where T<:AF

    # extrapolate the width up to the continuum
    intt2[end] = 1.0

    # interpolate bisector onto common intensity grid
    itp = linear_interp(intt1, bist, bc=last(bist))
    bist .= itp.(intt2)
    intt1 .= intt2

    # get indices to exclude data from fit
    thresh = 0.9
    idx1 = findfirst(x -> x .>= thresh, intt1)
    idx1 -= 1
    idx1 = clamp(idx1, firstindex(weights), lastindex(weights))

    depth = 1.0 - minimum(intt1)
    if depth > 0.45
        thresh = minimum(intt1) + 0.5 * depth
        idx2 = findfirst(x -> x .>= thresh, intt1)
        idx2 -= 1
        idx2 = clamp(idx2, firstindex(weights), idx1)
    else
        idx2 = 1
    end

    # set weights
    weights[1:2] .= 0.0
    weights[idx1:end] .= 0.0
    weights[1:idx2] .= 0.0

    # perform a polynomial fits to the bisector
    if (mu < 0.4) || (depth < 0.3)
        order = 1
    else
        order = 3
    end
    bfit1 = pfit(intt2, bist, order, weights=weights)
    bfit2 = pfit(intt2[3:10], bist[3:10], 1)

    # replace top and bottom with model fit
    bist[idx1:end] .= bfit1.(intt2[idx1:end])
    bist[1:3] .= bfit2.(intt2[1:3])
    return nothing
end

#=function extrapolate_input_data(bist::AA{T,1}, intt1::AA{T,1},
                                widt::AA{T,1}, intt2::AA{T,1},
                                mu::T; weights=ones(length(bist))) where T<:AF

    # extrapolate the width up to the continuum
    intt2[end] = 1.0

    # interpolate bisector onto common intensity grid
    itp = linear_interp(intt1, bist, bc=last(bist))
    bist .= itp.(intt2)
    intt1 .= intt2

    # get indices to exclude data from fit
    thresh = 0.9
    idx1 = findfirst(x -> x .>= thresh, intt1)
    idx1 -= 1
    idx1 = clamp(idx1, firstindex(weights), lastindex(weights))

    depth = 1.0 - minimum(intt1)
    if depth > 0.45
        thresh = minimum(intt1) + 0.5 * depth
        idx2 = findfirst(x -> x .>= thresh, intt1)
        idx2 -= 1
        idx2 = clamp(idx2, firstindex(weights), idx1)
    else
        idx2 = 1
    end

    # set weights
    weights[1:2] .= 0.0
    weights[idx1:end] .= 0.0
    weights[1:idx2] .= 0.0

    if (mu < 0.4) || (depth <= 0.3)
        order = 1
    else
        order = 2
    end
    bfit1 = pfit(intt2, bist, order, weights=weights)
    bfit2 = pfit(intt2[3:10], bist[3:10], 1)

    # replace top and bottom with model fit
    idx3 = findfirst(x -> x .>= 0.85, intt1)
    bist[idx3:end] .= bfit1.(intt2[idx3:end])
    bist[1:3] .= bfit2.(intt2[1:3])
    return nothing
end=#

function extrapolate_input_data(bis::AA{T,2}, int1::AA{T,2}, wid::AA{T,2}, int2::AA{T,2}, mu::T) where T<:AF
    weights = ones(size(bis,1))
    for t in 1:size(bis,2)
        # reset weights
        weights .= one(T)

        # take a slice for one time snapshot
        bist = view(bis, :, t)
        widt = view(wid, :, t)
        intt1 = view(int1, :, t)
        intt2 = view(int2, :, t)

        # plt.plot(bist, intt1)
        extrapolate_input_data(bist, intt1, widt, intt2, mu, weights=weights)
    end
    return nothing
end
