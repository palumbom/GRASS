abstract type AbstractSolarData end
struct SolarData{T1<:AF} <: AbstractSolarData
    wav::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    bis::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    dep::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    wid::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    len::Dict{Tuple{Symbol,Symbol}, Int}
end

"""
    SolarData(; dir=soldir, relative=true, ...)

Construct a `SpecParams` composite type instance.

# Arguments
- `dir::String=soldir`: Directory containing the pre-processed input data. Default directory is set in src/config.jl
- `relative::Bool=true`: Set whether wavelengths are on absolute scale or expressed relative to rest wavelength.
"""
function SolarData(;dir::String=soldir*"FeI_5434/", λrest::Float64=NaN,
                   relative::Bool=true, fixed_width::Bool=false,
                   fixed_bisector::Bool=false, extrapolate::Bool=true,
                   contiguous_only::Bool=false, adjust_mean::Bool=true)
    @assert isdir(dir)
    df = sort_input_data(dir=dir)
    return SolarData(df, λrest=λrest, relative=relative, fixed_width=fixed_width,
                     fixed_bisector=fixed_bisector, extrapolate=extrapolate,
                     contiguous_only=contiguous_only, adjust_mean=adjust_mean)
end

function SolarData(df::DataFrame; λrest::Float64=NaN,
                   relative::Bool=true, fixed_width::Bool=false,
                   fixed_bisector::Bool=false, extrapolate::Bool=true,
                   contiguous_only::Bool=false, adjust_mean::Bool=true)
    # allocate memory for data dictionaries
    wavdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    bisdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    depdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    widdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    lengths = Dict{Tuple{Symbol,Symbol}, Int}()

    # get rest wavelength for data if not passed
    if isnan(λrest)
        filename = glob("*_line_properties.h5", df.fpath[1])
        λrest = h5open(filename[1], "r") do f
            g = f["properties"]
            attr = HDF5.attributes(g)
            λrest = read(attr["air_wavelength"])
            return λrest
        end
    end

    # loop over unique mu + axis pairs
    for ax in unique(df.axis)
        for mu in unique(df.mu)
            # filter data for mu + axis pair
            f1 = x -> x == ax
            f2 = y -> y == mu
            df_temp = filter(:axis => f1, df)
            df_temp = filter(:mu => f2, df_temp)

            # move on if no data for mu + axis pair
            if isempty(df_temp)
                continue
            end

            # stitch the time series
            wavall, bisall, depall, widall = stitch_time_series(df_temp, adjust_mean=adjust_mean, contiguous_only=contiguous_only)

            # clean the input
            wavall, bisall, depall, widall = clean_input(wavall, bisall, depall, widall)

            # deal with NaNs in wavall measurements
            if extrapolate
                extrapolate_bisector(wavall, bisall)
            end

            # assign key-value pairs to dictionary
            if relative
                relwav = relative_bisector_wavelengths(wavall, λrest)
                wavdict[(Symbol(ax), Symbol(mu))] = relwav
            else
                wavdict[(Symbol(ax), Symbol(mu))] = wavall
            end
            bisdict[(Symbol(ax), Symbol(mu))] = bisall

            # assign width
            if fixed_width
                dep, wid = calc_fixed_width()
                depdict[(Symbol(ax), Symbol(mu))] = zeros(size(depall))
                widdict[(Symbol(ax), Symbol(mu))] = zeros(size(depall))
                depdict[(Symbol(ax), Symbol(mu))] .= dep
                widdict[(Symbol(ax), Symbol(mu))] .= wid
            else
                depdict[(Symbol(ax), Symbol(mu))] = depall
                widdict[(Symbol(ax), Symbol(mu))] = widall
            end

            # set bisector to zero if fixed
            if fixed_bisector
                wavdict[(Symbol(ax), Symbol(mu))] = zeros(size(wavall))
            end

            # get length
            lengths[(Symbol(ax), Symbol(mu))] = size(wavall,2)
        end
    end
    return SolarData(wavdict, bisdict, depdict, widdict, lengths)
end

function calc_fixed_width(;λ::T1=5434.5232, M::T1=26.0, T::T1=5778.0, v_turb::T1=3.5e5) where T1<:AF
    wid = width_thermal(λ=λ, M=M, T=T, v_turb=v_turb)
    λs = collect(range(λ - 1.0, λ + 1.0, step=λ/7e5))
    prof = gaussian_line.(λs, mid=λ, width=wid, depth=0.86)
    depall, widall = calc_width_at_depth(λs, prof, center=λ, len=100)
    widall[1] = 0.01
    return depall, widall
end

function clean_input(wavall::AA{T,2}, bisall::AA{T,2}, depall::AA{T,2}, widall::AA{T,2}) where T<:AF
    @assert size(wavall) == size(bisall)
    @assert size(depall) == size(widall)

    # make boolean array (column will be stripped if badcol[i] == true)
    badcols = zeros(Bool, size(wavall,2))

    # find spread of data
    wav_std = std(wavall, dims=2)
    wid_std = std(wavall, dims=2)

    # loop through checking for bad columns
    for i in 1:size(wavall,2)
        # check for monotinicity
        if !ismonotonic(widall[:,i])
            badcols[i] = true
        end

        # check for bad width measurements
        if all(iszero(widall[:,i]))
            badcols[i] = true
        end

        # check for excessive NaNs
        idx = findfirst(x->isnan(x), wavall[:, i])
        if bisall[idx,i] < 0.85
            badcols[i] = true
        end

        # remove data that is significant outlier (bisector)
        if any(abs.(wavall[5:50,1] - wavall[5:50,i]) .> (4.0 .* wav_std[5:50]))
            badcols[i] = true
        end
    end

    # strip the bad columns and return new arrays
    return strip_columns(wavall, badcols), strip_columns(bisall, badcols),
           strip_columns(depall, badcols), strip_columns(widall, badcols)
end

function extrapolate_bisector(wavall::AA{T,2}, bisall::AA{T,2}; top::T=0.8) where T<:AF
    for i in 1:size(wavall,2)
        # take a slice for one time snapshot
        wavt = view(wavall, :, i)
        bist = view(bisall, :, i)

        # find the last index for good data
        ind1 = searchsortedfirst(bist, top)
        ind2 = searchsortedfirst(bist, top-0.1)

        # calculate an average slope and project it up to continuum
        dydx = mean((wavt[(ind2+1:ind1)] .- wavt[(ind2:ind1-1)]) ./ (bist[ind2+1:ind1] .- bist[(ind2:ind1-1)]))
        wavall[ind1:end, i] .= (dydx .* (bist[ind1:end] .- bist[ind1]) .+ wavt[ind1])
    end
    return nothing
end

function relative_bisector_wavelengths(wav::AA{T,2}, λrest::T) where T<:AF
    # allocate array of size(wav)
    relwav = similar(wav)
    for i in 1:size(wav,2)
        λgrav = (635.0/c_ms) * λrest
        relwav[:,i] = wav[:,i] .- λrest .- λgrav
    end
    return relwav
end
