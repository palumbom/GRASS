struct SolarData{T1<:AF}
    bis::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    int::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    wid::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    top::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,1}}
    dep_contrast::Dict{Tuple{Symbol,Symbol}, T1}
    cbs::Dict{Tuple{Symbol,Symbol}, T1}
    len::Dict{Tuple{Symbol,Symbol}, Int}
    ax::Array{Symbol,1}
    mu::Array{Symbol,1}
end

"""
    SolarData(; dir=soldir, relative=true, ...)

Construct a `SpecParams` composite type instance.

# Arguments
- `dir::String=soldir`: Directory containing the pre-processed input data. Default directory is set in src/config.jl
- `relative::Bool=true`: Set whether bisectors are measured on absolute wavelength scale or expressed relative to rest wavelength of the line.
"""
function SolarData(;fname::String="", relative::Bool=true, extrapolate::Bool=true,
                    adjust_mean::Bool=true, contiguous_only::Bool=false,
                    fixed_width::Bool=false, fixed_bisector::Bool=false)
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
                     fixed_width=fixed_width, fixed_bisector=fixed_bisector)
end

function SolarData(fname::String; relative::Bool=true, extrapolate::Bool=true,
                   adjust_mean::Bool=true, contiguous_only::Bool=false,
                   fixed_width::Bool=false, fixed_bisector::Bool=false)
    # make sure the file exists
    @assert isfile(fname)

    # initialize data structure for input data
    topdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,1}}()
    bisdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    intdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    widdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    cbsdict = Dict{Tuple{Symbol,Symbol}, Float64}()
    lendict = Dict{Tuple{Symbol,Symbol}, Int}()
    axs = []
    mus = []

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
            badcols = identify_bad_cols(bis, int1, wid, int2)
            if sum(badcols) > 0
                top = strip_columns(top, badcols)
                bis = strip_columns(bis, badcols)
                wid = strip_columns(wid, badcols)
                int1 = strip_columns(int1, badcols)
                int2 = strip_columns(int2, badcols)
            end

            if extrapolate
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
            push!(axs, ax)
            push!(mus, mu)
        end
    end

    # get key components
    axs = Symbol.(unique(axs))
    mus = Symbol.(sort!(unique(mus)))

    # compute mean depth at disk center
    dep_dc = mean(1.0 .- view(intdict[(:c, :mu10)], 1, :))

    # compute depth contrasts
    dep_contrast = Dict{Tuple{Symbol,Symbol}, Float64}()
    for k in keys(intdict)
        dep_contrast[k] = mean(1.0 .- view(intdict[k], 1, :)) / dep_dc
    end

    # make sure the keys of all dictionaries are sorted
    bisdict = Dict(sort(bisdict, by=x -> x[2]))
    intdict = Dict(sort(intdict, by=x -> x[2]))
    widdict = Dict(sort(widdict, by=x -> x[2]))
    topdict = Dict(sort(topdict, by=x -> x[2]))
    dep_contrast = Dict(sort(dep_contrast, by=x -> x[2]))
    cbsdict = Dict(sort(cbsdict, by=x -> x[2]))
    lendict = Dict(sort(lendict, by=x -> x[2]))

    # assertion to verify keys are the same
    @assert all(keys(bisdict) .== keys(dep_contrast))

    # construct the composite tpye
    return SolarData(bisdict, intdict, widdict, topdict, dep_contrast, cbsdict, lendict, axs, mus)
end
