struct SpecParams{T<:AF}
    lines::AA{T,1}
    depths::AA{T,1}
    conv_blueshifts::AA{T,1}
    variability::AA{Bool,1}
    coverage::Tuple{T,T}
    resolution::T
    lambdas::AA{T,1}
    soldata::SolarData{T}
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
function SpecParams(;lines=[], depths=[], variability=[], resolution=7e5,
                    buffer=0.75, fixed_width=false, fixed_bisector=false,
                    extrapolate=true, contiguous_only=false)
    @assert !isempty(lines)
    @assert !isempty(depths)
    @assert length(lines) == length(depths)
    @assert buffer >= 0.75

    # assign depths if needed, otherwise check same number of centers + depths
    if any(isnan.(depths))
        depths = rand(length(lines))
    else
        @assert !any(depths .> 1.0)
        @assert !any(depths .< 0.0)
    end

    # read in convective blueshifts
    df = CSV.read(soldir * "../convective_blueshift.dat", DataFrame,
                  header=1, delim=",", type=Float64)

    # assign convective blueshifts
    blueshifts = similar(lines)
    for i in eachindex(depths)
        ind = searchsortednearest(df.depth, depths[i])
        blueshifts[i] = df.blueshift[ind]
    end

    # assign fixed_width booleans
    if isempty(variability)
        variability = trues(length(lines))
    else
        @assert length(lines) == length(variability)
    end

    # calculate coverage
    # TODO: constant delta v, not wavelength?
    coverage = (minimum(lines) - buffer, maximum(lines) + buffer)

    # get generator for lambdas
    dλ = coverage[1]/resolution
    lambdas = coverage[1]:dλ:coverage[2]

    # get observational data
    soldata = SolarData(relative=true, fixed_width=fixed_width,
                        fixed_bisector=fixed_bisector, extrapolate=extrapolate,
                        contiguous_only=contiguous_only)
    return SpecParams(lines, depths, blueshifts, variability, coverage, resolution, lambdas, soldata)
end
