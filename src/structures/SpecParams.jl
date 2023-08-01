# get groups of lines in same osbserved spectral regions
const line_groups = [["FeI_5250.2", "FeI_5250.6"],
                     ["FeI_5379", "CI_5380", "TiII_5381", "FeI_5382", "FeI_5383"],
                     ["MnI_5432", "FeI_5432", "FeI_5434", "NiI_5435", "FeI_5436.3", "FeI_5436.6"],
                     ["FeI_5576", "NiI_5578"],
                     ["NaI_5896"],
                     ["FeII_6149", "FeI_6151"],
                     ["CaI_6169.0", "CaI_6169.5", "FeI_6170", "FeI_6173"],
                     ["FeI_6301", "FeI_6302"]]

function get_name_from_filename(line1::String)
    # get filename from full path
    split1 = splitdir(line1)[end]

    # get name if line is filenames
    split1 = split(split1, ".h5")
    return split1[1]
end

function in_same_group(line1::String, line2::String)
    # get name if lines are filenames
    split1 = get_name_from_filename(line1)
    split2 = get_name_from_filename(line2)

    # check if they are in the same group
    for row in line_groups
        if split1 in row
            return split2 in row
        else
            continue
        end
    end
    return nothing
end

function get_template_wavelength(line1::String)
    # get filename
    if split(line1, ".")[end] != "h5"
        fname = line1 * ".h5"
    else
        fname = line1
    end

    # open the file and get the rest wavelength
    λrest = h5open(fname, "r") do f
        # get rest wavelength for line
        attr = HDF5.attributes(f)
        λrest = read(attr["air_wavelength"])
    end
    return λrest
end

"""
    SpecParams(lines, depths, geffs, conv_blueshifts, variability, resolution, lambdas, templates)
"""
struct SpecParams{T<:AF}
    lines::AA{T,1}
    depths::AA{T,1}
    geffs::AA{T,1}
    conv_blueshifts::AA{T,1}
    variability::AA{Bool,1}
    templates::AA{String,1}
    resolution::T
    lambdas::AA{T,1}
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
    @assert buffer > 0.0

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
            # blueshifts[i] = rand(Normal(df.blueshift[idx], df.sigma[idx]))
            blueshifts[i] = rand(Normal(df.blueshift[idx], 0.0))
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
        @warn "No line template specified!"

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

    # get indices to sort on template line wavelength
    template_wavelengths = get_template_wavelength.(templates)
    inds = sortperm(template_wavelengths)

    # collect ranges if necessary
    lines = collect(lines)
    depths = collect(depths)
    geffs = collect(geffs)
    blueshifts = collect(blueshifts)

    # now do the sorting and return
    return SpecParams(lines[inds], depths[inds], geffs[inds], blueshifts[inds],
                      variability[inds], templates[inds], resolution, lambdas)
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
        file = joinpath(GRASS.soldir, file)
    end
    @assert isfile(file)

    # get indices
    idx = spec.templates .== file
    return SpecParams(view(spec.lines, idx),
                      view(spec.depths, idx),
                      view(spec.geffs, idx),
                      view(spec.conv_blueshifts, idx),
                      view(spec.variability, idx),
                      view(spec.templates, idx),
                      spec.resolution,
                      spec.lambdas)
end
