"""
    SolarData

Container for line profile data and metadata loaded from a template HDF5 file.
"""
struct SolarData{T1<:AF}
    bis::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    int::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    wid::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    top::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,1}}
    dep_contrast::Dict{Tuple{Symbol,Symbol}, T1}
    cbs::Dict{Tuple{Symbol,Symbol}, T1}
    len::Dict{Tuple{Symbol,Symbol}, Int}
    ax::Array{Int,1}
    mu::Array{T1,1}
end

"""
    SolarData(; fname="", relative=true, extrapolate=true, adjust_mean=true,
              contiguous_only=false, fixed_width=false, fixed_bisector=false,
              strip_cols=true)

Load solar line profile data from a template HDF5 file.

# Keyword Arguments
- `fname::String=""`: template filename; defaults to `soldir/FeI_5434.h5`.
- `relative::Bool=true`: store bisectors relative to the rest wavelength.
- `extrapolate::Bool=true`: extrapolate bisectors when intensities fall below 0.75.
- `adjust_mean::Bool=true`: match means across noncontiguous datasets.
- `contiguous_only::Bool=false`: load only the first contiguous time series per disk position.
- `fixed_width::Bool=false`: use a fixed width for all epochs at a disk position.
- `fixed_bisector::Bool=false`: use a fixed bisector for all epochs at a disk position.
- `strip_cols::Bool=true`: remove columns flagged as bad before processing.
"""
function SolarData(;fname::String="", relative::Bool=true, extrapolate::Bool=true,
                    adjust_mean::Bool=true, contiguous_only::Bool=false,
                    fixed_width::Bool=false, fixed_bisector::Bool=false,
                    strip_cols::Bool=true)
    # default
    if isempty(fname)
        fname = joinpath(soldir, "FeI_5434.h5")
    end

    # if doesn't end in file
    if !contains(fname, ".h5")
        fname *= ".h5"
    end

    if !isabspath(fname)
        fname = joinpath(soldir, fname)
    end

    return SolarData(fname, relative=relative, extrapolate=extrapolate,
                     adjust_mean=adjust_mean, contiguous_only=contiguous_only,
                     fixed_width=fixed_width, fixed_bisector=fixed_bisector,
                     strip_cols=strip_cols)
end

function SolarData(fname::String; relative::Bool=true, extrapolate::Bool=true,
                   adjust_mean::Bool=true, contiguous_only::Bool=false,
                   fixed_width::Bool=false, fixed_bisector::Bool=false,
                   strip_cols::Bool=true)
    # make sure the file exists
    @assert isfile(fname)

    # initialize data structure for input data
    topdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,1}}()
    bisdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    intdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    widdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    cbsdict = Dict{Tuple{Symbol,Symbol}, Float64}()
    lendict = Dict{Tuple{Symbol,Symbol}, Int}()

    # open the file
    h5open(fname, "r") do f
        # get rest wavelength for line
        attr = HDF5.attributes(f)
        λrest = read(attr["air_wavelength"])

        # loop over the keys corresponding to disk positions
        for k in keys(f)
            # get the axis and mu values
            attr = HDF5.attributes(f[k])
            ax = read(attr["axis"])
            mu = read(attr["mu"])
            vconv = read(attr["vconv"]) / c_ms # convert to z = v/c

            # only read in the first set of observations
            if contiguous_only
                # get key and length of data set
                t = first(keys(f[k]))
                ntimes = [read(HDF5.attributes(f[k][t])["length"])]

                # read in
                top = read(f[k][t]["top_ints"])
                bis = read(f[k][t]["bisectors"])
                int1 = read(f[k][t]["bis_intensities"])
                wid = read(f[k][t]["widths"])
                int2 = read(f[k][t]["wid_intensities"])
            # stitch together all observations of given disk position
            else
                # get total number of epochs
                ntimes = zeros(Int, length(f[k]))
                for (i, t) in enumerate(keys(f[k]))
                    ntimes[i] = read(HDF5.attributes(f[k][t])["length"])
                end

                # allocate memory
                top = zeros(sum(ntimes))
                bis = zeros(100, sum(ntimes))
                wid = zeros(100, sum(ntimes))
                int1 = zeros(100, sum(ntimes))
                int2 = zeros(100, sum(ntimes))

                # read in
                for (i, t) in enumerate(keys(f[k]))
                    top[sum(ntimes[1:i-1])+1:sum(ntimes[1:i])] = read(f[k][t]["top_ints"])
                    bis[:, sum(ntimes[1:i-1])+1:sum(ntimes[1:i])] = read(f[k][t]["bisectors"])
                    int1[:, sum(ntimes[1:i-1])+1:sum(ntimes[1:i])] = read(f[k][t]["bis_intensities"])
                    wid[:, sum(ntimes[1:i-1])+1:sum(ntimes[1:i])] = read(f[k][t]["widths"])
                    int2[:, sum(ntimes[1:i-1])+1:sum(ntimes[1:i])] = read(f[k][t]["wid_intensities"])
                end
            end

            # identify bad columns and strip them out
            if strip_cols
                badcols = identify_bad_cols(bis, int1, wid, int2)
                if sum(badcols) > 0
                    top = strip_columns(top, badcols)
                    bis = strip_columns(bis, badcols)
                    wid = strip_columns(wid, badcols)
                    int1 = strip_columns(int1, badcols)
                    int2 = strip_columns(int2, badcols)
                end
            end

            if extrapolate & (minimum(int1) < 0.75)# & !(contains(fname, "FeI_5383"))
                extrapolate_input_data(bis, int1, wid, int2, parse_mu_string(mu))
            end

            # match the means of noncontiguous datasets
            if length(ntimes) > 1 && adjust_mean && !contiguous_only
                # fix ntimes to deal with removed columns
                inds = findall(badcols)
                ntimes_cum = cumsum(ntimes)
                for i in inds
                    idx = findfirst(x -> x .>= i, ntimes_cum)
                    ntimes[idx] -= 1
                end

                # now adjust the mean
                adjust_data_mean(bis, ntimes)
            end

            # express bisector wavelengths relative to mean
            # add in convective blueshift at line synthesis stage
            if relative
                relative_bisector_wavelengths(bis)
            end

            # assign input data to dictionary
            lendict[Symbol(ax), Symbol(mu)] = size(bis, 2)
            topdict[Symbol(ax), Symbol(mu)] = top
            cbsdict[Symbol(ax), Symbol(mu)] = vconv
            bisdict[Symbol(ax), Symbol(mu)] = (bis .* !fixed_bisector) .+ (λrest .* fixed_bisector * !relative)
            intdict[Symbol(ax), Symbol(mu)] = int1
            if !fixed_width
                widdict[Symbol(ax), Symbol(mu)] = wid
            else
                wid .= wid[:,1]
                widdict[Symbol(ax), Symbol(mu)] = wid
            end
        end
    end

    # compute mean depth at disk center
    dep_dc = mean(1.0 .- view(intdict[(:c, :mu10)], 1, :))

    # compute depth contrasts
    dep_contrast = Dict{Tuple{Symbol,Symbol}, Float64}()
    for k in keys(intdict)
        dep_contrast[k] = mean(1.0 .- view(intdict[k], 1, :)) / dep_dc
    end

    # make sure the keys of all dictionaries are sorted the same way
    bisdict = Dict(sort!(OrderedDict(bisdict), by=x -> x[2]))
    intdict = Dict(sort!(OrderedDict(intdict), by=x -> x[2]))
    widdict = Dict(sort!(OrderedDict(widdict), by=x -> x[2]))
    topdict = Dict(sort!(OrderedDict(topdict), by=x -> x[2]))
    dep_contrast = Dict(sort!(OrderedDict(dep_contrast), by=x -> x[2]))
    cbsdict = Dict(sort!(OrderedDict(cbsdict), by=x -> x[2]))
    lendict = Dict(sort!(OrderedDict(lendict), by=x -> x[2]))

    # assertion to verify keys are the same
    @assert all(keys(bisdict) .== keys(dep_contrast))

    # get mu and ax codes
    disc_ax = parse_ax_string.(getindex.(keys(lendict),1))
    disc_mu = parse_mu_string.(getindex.(keys(lendict),2))

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

    # construct the composite tpye
    return SolarData(bisdict, intdict, widdict, topdict, dep_contrast, cbsdict, lendict, disc_ax, disc_mu)
end
