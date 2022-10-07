struct SpecParams{T<:AF}
    lines::AA{T,1}
    depths::AA{T,1}
    geffs::AA{T,1}
    conv_blueshifts::AA{T,1}
    variability::AA{Bool,1}
    resolution::T
    lambdas::AA{T,1}
    templates::AA{String,1}
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
                    templates=[], blueshifts=[], resolution=7e5, buffer=2.0)
    @assert length(lines) == length(depths)
    @assert !isempty(lines)
    @assert !isempty(depths)
    @assert all(depths .< 1.0) && all(depths .> 0.0)
    @assert buffer >= 1.0

    # assign lande g factors if they haven't been
    if isempty(geffs)
        geffs = zeros(length(lines))
    end

    # read in convective blueshifts
    df = CSV.read(datdir * "convective_blueshift.dat", DataFrame,
                  header=1, delim=",", types=Float64)

    # assign convective blueshifts and
    if isempty(blueshifts)
        blueshifts = similar(lines)
        for i in eachindex(depths)
            idx = searchsortednearest(df.depth, depths[i])
            blueshifts[i] = df.blueshift[idx]
        end
    else
        @assert length(blueshifts) == length(lines)
   end

   # convert blueshift to z=v/c
   blueshifts ./=  c_ms

    # assign fixed_width booleans
    if isempty(variability)
        variability = trues(length(lines))
    else
        @assert length(lines) == length(variability)
    end

    # generate Delta ln lambda
    minλ = minimum(lines) - buffer
    maxλ = maximum(lines) + buffer
    Δlnλ = (1.0 / resolution)
    lnλs = range(log(minλ), log(maxλ), step=Δlnλ)
    lambdas = exp.(lnλs)

    # assign each line to the input data to synth it from
    # TODO move this to its own function and make it more thought out
    if isempty(templates)
        # get properties of input data lines
        lp = LineProperties()
        geff_input = get_geff(lp)
        depth_input = get_depth(lp)
        input_files = get_file(lp)

        # allocate memory and loop, choosing best template line
        templates = Array{String, 1}(undef, length(lines))
        for i in eachindex(templates)
            param_dist = sqrt.((geff_input .- geffs[i]).^2 + (depth_input .- depths[i]).^2)
            templates[i] = input_files[argmin(param_dist)]
        end
    else
        @assert length(templates) == length(lines)
        if all(map(x -> split(x, ".")[end] != "h5", templates))
            templates .*= ".h5"
        end
    end

    # make sure templates are absolute paths to files
    if all(.!isabspath.(templates))
        for i in eachindex(templates)
            templates[i] = GRASS.soldir * templates[i]
        end
    end
    @assert all(isfile.(templates))

    # now make sure everything is sorted
    if !issorted(lines)
        inds = sortperm(lines)
        lines = view(lines, inds)
        depths = view(depths, inds)
        geffs = view(geffs, inds)
        blueshifts = view(blueshifts, inds)
        variability = view(variability, inds)
        templates = view(templates, inds)
    end
    return SpecParams(lines, depths, geffs, blueshifts,
                      variability, resolution,
                      lambdas, templates)
end

"""
    SpecParams(spec, template_file)

Return a copy of ```spec``` with only the synthetic line parameters corresponding to ```template_file```
"""
function SpecParams(spec::SpecParams, template_file::String)
    # make sure it's a file name and not just a string
    file = template_file
    if split(file, ".")[end] != "h5"
        file *= ".h5"
    end

    # make sure it's an absolute path
    if !isabspath(file)
        file = GRASS.soldir * file
    end
    @assert isfile(file)

    # get indices
    idx = findall(spec.templates .== file)
    return SpecParams(view(spec.lines, idx),
                      view(spec.depths, idx),
                      view(spec.geffs, idx),
                      view(spec.conv_blueshifts, idx),
                      view(spec.variability, idx),
                      spec.resolution,
                      view(spec.lambdas, :),
                      view(spec.templates, idx))
end
