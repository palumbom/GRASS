# abstract type AbstractSpecParams end
struct SpecParams{T<:AF} #<: AbstractSpecParams
    lines::AA{T,1}
    depths::AA{T,1}
    geffs::AA{T,1}
    conv_blueshifts::AA{T,1}
    variability::AA{Bool,1}
    coverage::Tuple{T,T}
    resolution::T
    lambdas::AA{T,1}
    indata::InputData
    data_inds::AA{Int64,1}
    kwargs::Base.Pairs
end


"""
    SpecParams(; lines=[], depths=[], variability=[], resolution=7e5)

Construct a `SpecParams` composite type instance. If `variability` is not specified, all lines are variable by default.

# Arguments
- `lines::AbstractArray{Float64,1}=[]`: List of line centers (in angstroms)
- `depths::AbstractArray{Float64,1}=[]`: List of line depths
- `variability::AbstractArray{Bool,1}=[]`: Array of booleans controlling whether corresponding line has variability.
- `resolution::Float64=7e8`: Spectral resolution of spectrum
"""
function SpecParams(;lines=[], depths=[], geffs=[], variability=[],
                    indirs=[], resolution=7e5, buffer=1.0, kwargs...)
    @assert length(lines) == length(depths)
    @assert !isempty(lines)
    @assert !isempty(depths)
    @assert all(depths .< 1.0)
    @assert all(depths .> 0.0)
    @assert buffer >= 0.75

    # assign lande g factors if they haven't been
    if isempty(geffs)
        geffs = zeros(length(lines))
    end

    # read in convective blueshifts
    df = CSV.read(soldir * "/convective_blueshift.dat", DataFrame,
                  header=1, delim=",", type=Float64)

    # assign convective blueshifts and convert blueshift to z=v/c
    blueshifts = similar(lines)
    for i in eachindex(depths)
        idx = searchsortednearest(df.depth, depths[i])
        blueshifts[i] = df.blueshift[idx] / c_ms
    end

    # assign fixed_width booleans
    if isempty(variability)
        variability = trues(length(lines))
    else
        @assert length(lines) == length(variability)
    end

    # calculate coverage
    coverage = (minimum(lines) - buffer, maximum(lines) + buffer)

    # generate Delta ln lambda
    Δlnλ = (1.0 / resolution)
    lnλs = range(log(coverage[1]), log(coverage[2]), step=Δlnλ)
    lambdas = exp.(lnλs)

    # tabulate all available input data
    indata = InputData()

    # assign indices that point synthetic line to appropriate input data
    geff_input = get_geff.(indata.lineprops)
    depth_input = get_depth.(indata.lineprops)
    data_inds = zeros(Int64, length(lines))

    # loop over lines and do 2D nearest neighbor
    if isempty(indirs)
        for i in eachindex(lines)
            param_dist = sqrt.((geff_input .- geffs[i]).^2 + (depth_input .- depths[i]).^2)
            idx = argmin(param_dist)
            data_inds[i] = idx
        end
    else
        @assert length(indirs) == length(lines)
        for i in eachindex(indirs)
            data_inds[i] = findfirst(indirs[i] .== indata.dirs)
        end
    end

    # now make sure everything is sorted
    if !issorted(lines)
        inds = sortperm(lines)
        lines = lines[inds]
        depths = depths[inds]
        geffs = geffs[inds]
        blueshifts = blueshifts[inds]
        variability = variability[inds]
        data_inds = data_inds[inds]
    end
    return SpecParams(lines, depths, geffs, blueshifts,
                      variability, coverage, resolution,
                      lambdas, indata, data_inds, kwargs)
end

function SpecParams(spec::SpecParams, idx::Int64)
    inds = spec.data_inds .== idx
    indata_temp = InputData(spec.indata.dirs[idx], spec.indata.lineprops[idx])
    return SpecParams(spec.lines[inds], spec.depths[inds], spec.geffs[inds],
                      spec.conv_blueshifts[inds], spec.variability[inds],
                      spec.coverage, spec.resolution, spec.lambdas,
                      indata_temp, spec.data_inds[inds], spec.kwargs)
end

function SpecParams(config::String)
    @assert isfile(config)



    return nothing
end
