abstract type AbstractSolarData end
struct SolarData{T1<:AF} <: AbstractSolarData
    wav::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    bis::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    dep::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    wid::Dict{Tuple{Symbol,Symbol}, AbstractArray{T1,2}}
    len::Dict{Tuple{Symbol,Symbol}, Int}
    ax::Array{Symbol,1}
    mu::Array{Symbol,1}
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

function SolarData(fname::String, relative::Bool=true, fixed_width::Bool=false,
                   fixed_bisector::Bool=false, extrapolate::Bool=true,
                   contiguous_only::Bool=false, adjust_mean::Bool=true)
    # set up data structure for input data
    wavdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    bisdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    depdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    widdict = Dict{Tuple{Symbol,Symbol}, AA{Float64,2}}()
    lengths = Dict{Tuple{Symbol,Symbol}, Int}()

    # get rest wavelength for line
    λrest = h5open(fname, "r") do f
        attr = HDF5.attributes(f)
        λrest = read(attr["air_wavelength"])
        return λrest
    end

    # loop over unique mu + axis pairs
    axs = unique(df.axis)
    mus = unique(df.mu)
    for ax in axs
        for mu in mus
            # filter the data set
            f3 = (x,y) -> ((x == ax) & (y == mu))
            df_temp = filter([:axis, :mu] => f3, df)

            # move on if no data for mu + axis pair
            if isempty(df_temp)
                continue
            end

            # stitch the time series
            wavall, bisall, depall, widall = stitch_time_series(df_temp, adjust_mean=adjust_mean, contiguous_only=contiguous_only)

            # clean the input (loop through twice to catch outliers)
            for i in 1:2
                wavall, bisall, depall, widall = clean_input(wavall, bisall, depall, widall)
            end

            # extrapolate over data where uncertainty explodes
            if extrapolate
                if tryparse(Int64, mu[end-1:end]) <= 6 # ie mu <= 0.5
                    top = 0.8
                else
                    top = 0.9
                end
                extrapolate_bisector(wavall, bisall, top=top)
                extrapolate_width(depall, widall)
            end

            # assign key-value pairs to dictionary
            if relative
                relative_bisector_wavelengths(wavall, λrest)
            end
            wavdict[(Symbol(ax), Symbol(mu))] = wavall
            bisdict[(Symbol(ax), Symbol(mu))] = bisall

            # assign width
            if fixed_width
                dep, wid = calc_fixed_width()
                depdict[(Symbol(ax), Symbol(mu))] = zeros(size(depall))
                widdict[(Symbol(ax), Symbol(mu))] = zeros(size(depall))
                depdict[(Symbol(ax), Symbol(mu))] .= dep
                widdict[(Symbol(ax), Symbol(mu))] .= wid

                extrapolate_width(depdict[(Symbol(ax), Symbol(mu))],
                                  widdict[(Symbol(ax), Symbol(mu))])
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
    return SolarData(wavdict, bisdict, depdict, widdict, lengths, Symbol.(axs), sort!(Symbol.(mus)))
end

function clean_input(wavall::AA{T,2}, bisall::AA{T,2}, depall::AA{T,2}, widall::AA{T,2}) where T<:AF
    @assert size(wavall) == size(bisall)
    @assert size(depall) == size(widall)

    # make boolean array (column will be stripped if badcol[i] == true)
    badcols = zeros(Bool, size(wavall,2))

    # find spread of data
    wav_std = std(wavall, dims=2)
    wid_std = std(widall, dims=2)

    # find mean and median of data
    wav_avg = mean(wavall, dims=2)
    wid_avg = mean(widall, dims=2)
    wav_med = median(wavall, dims=2)
    wid_med = median(widall, dims=2)

    # loop through checking for bad columns
    for i in 1:size(wavall,2)
        wavt = view(wavall, :, i)
        bist = view(bisall, :, i)
        dept = view(depall, :, i)
        widt = view(widall, :, i)

        # check for monotinicity
        if !ismonotonic(widt[.!isnan.(widt)])
            badcols[i] = true
        end

        # check for bad width measurements
        if all(iszero(widt))
            badcols[i] = true
        end

        # check for excessive NaNs
        idx = findfirst(x->isnan(x), wavt)
        if !isnothing(idx) && (bist[idx] < 0.85)
            badcols[i] = true
        end

        # remove data that is significant outlier (bisector)
        wav_cond = any(abs.(wav_avg[5:50] .- wavt[5:50]) .> (4.0 .* wav_std[5:50]))
        wid_cond = any(abs.(wid_avg .- widt) .> (4.0 .* wid_std))
        if wav_cond || wid_cond
            badcols[i] = true
        end
    end

    # strip the bad columns and return new arrays
    return strip_columns(wavall, badcols), strip_columns(bisall, badcols),
           strip_columns(depall, badcols), strip_columns(widall, badcols)
end

function extrapolate_bisector(wavall::AA{T,2}, bisall::AA{T,2}; top::T=0.9) where T<:AF
    for i in 1:size(wavall,2)
        # take a slice for one time snapshot
        wavt = view(wavall, :, i)
        bist = view(bisall, :, i)

        # first fix the bottom-most measurements in bisector
        dydx = (wavt[5] - wavt[4])/(bist[5] - bist[4])
        wavt[1:3] .= dydx * (bist[1:3] .- bist[4]) .+ wavt[4]

        # find the last index for good data
        ind1 = searchsortedfirst(bist, top)
        ind2 = searchsortedfirst(bist, top-0.1)

        # calculate an average slope and extrapolate it up to continuum
        dydx = mean((wavt[(ind2+1:ind1)] .- wavt[(ind2:ind1-1)]) ./ (bist[ind2+1:ind1] .- bist[(ind2:ind1-1)]))
        wavt[ind2:end] .= (dydx .* (bist[ind2:end] .- bist[ind2]) .+ wavt[ind2])
    end
    return nothing
end

function extrapolate_width(depall::AA{T,2}, widall::AA{T,2}) where T<:AF
    for i in 1:size(depall,2)
        # take a slice for one time snapshot
        dept = view(depall, :, i)
        widt = view(widall, :, i)

        if !isnan(lastindex(widt))
            continue
        end

        # mask nans
        dept_mask = dept[.!isnan.(widt)]
        widt_mask = widt[.!isnan.(widt)]

        # fit a polynomial and extrapolate width up to continuum
        idx = findlast(x -> .!isnan.(x), widt)
        poly = pfit(dept_mask[idx-2:idx], widt_mask[idx-2:idx], 1)

        dept[end] = 1.0
        widt[idx:end] .= poly.(dept[idx:end])
    end
    return nothing
end

function relative_bisector_wavelengths(wav::AA{T,2}, λrest::T) where T<:AF
    λgrav = (635.0/c_ms) * λrest
    for i in 1:size(wav,2)
        wav[:,i] .-= (λrest .+ λgrav)
    end
    return nothing
end

function calc_fixed_width(;λ::T1=5434.5232, M::T1=26.0, T::T1=5778.0, v_turb::T1=3.5e5) where T1<:AF
    wid = width_thermal(λ=λ, M=M, T=T, v_turb=v_turb)
    λs = collect(range(λ - 1.0, λ + 1.0, step=λ/7e5))
    prof = gaussian_line.(λs, mid=λ, width=wid, depth=0.86)
    depall, widall = calc_width_function(λs, prof, nflux=100)
    widall[1] = 0.01
    return depall, widall
end
