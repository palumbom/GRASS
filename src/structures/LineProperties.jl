struct LineProperties{T<:AF}
    species::AA{String,1}
    mass::AA{T,1}             # amu
    depth::AA{T,1}
    λrest::AA{T,1}            # angstroms
    geff::AA{T,1}
    height::AA{T,1}           # km
    lower_level::AA{T,1}      # eV
    upper_level::AA{T,1}      # eV
    file::AA{String,1}
end

function LineProperties(;dir=GRASS.soldir, exclude::AA{String,1}=["CI_5380", "FeI_5382", "NaI_5896"])
    @assert isdir(dir)

    # get list of files
    files = glob("*.h5", dir)

    # allocate memory
    species = Array{String,1}()
    mass = Array{Float64,1}()
    depth = Array{Float64,1}()
    λrest = Array{Float64,1}()
    geff = Array{Float64,1}()
    height = Array{Float64,1}()
    lower_level = Array{Float64,1}()
    upper_level = Array{Float64,1}()
    filename = Array{String,1}()

    # loop over files
    for file in files
        # skip iteration if we want to skip this input data
        if splitdir(file)[end] in exclude
            println(">>> Excluding " * splitdir(file)[end])
            continue
        end

        # open the file and copy properties to arrays
        push!(filename, file)
        h5open(file, "r") do f
            # parse out the attributes
            attr = HDF5.attributes(f)
            push!(species, read(attr["species"]))
            push!(mass, read(attr["mass"]))
            push!(depth, read(attr["depth"]))
            push!(λrest, read(attr["air_wavelength"]))
            push!(geff, read(attr["g_eff"]))
            push!(height, read(attr["height"]))
            push!(lower_level, read(attr["lower_level"]))
            push!(upper_level, read(attr["upper_level"]))
        end
    end
    return LineProperties(species, mass, depth, λrest, geff, height, lower_level, upper_level, filename)
end

get_species(lp::LineProperties) = lp.species
get_rest_wavelength(lp::LineProperties) = lp.λrest
get_energy(lp::LineProperties) = lp.upper_level - lp.lower_level
get_upper_level(lp::LineProperties) = lp.upper_level
get_lower_level(lp::LineProperties) = lp.lower_level
get_depth(lp::LineProperties) = lp.depth
get_geff(lp::LineProperties) = lp.geff
get_height(lp::LineProperties) = lp.height
get_file(lp::LineProperties) = lp.file
get_name(lp::LineProperties) = map(x -> split(splitdir(x)[end], ".h5")[1], lp.file)
