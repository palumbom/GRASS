abstract type AbstractLineProperties end
struct LineProperties{T<:AF} <: AbstractLineProperties
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

function LineProperties(;dir=GRASS.soldir)
    @assert isdir(dir)

    # get list of files
    files = glob("*.h5", dir)

    # allocate memory
    species = Array{String,1}(undef, length(files))
    mass = zeros(length(files))
    depth = zeros(length(files))
    λrest = zeros(length(files))
    geff = zeros(length(files))
    height = zeros(length(files))
    lower_level = zeros(length(files))
    upper_level = zeros(length(files))
    filename = Array{String,1}(undef, length(files))

    # loop over files
    for i in eachindex(files)
        filename[i] = files[i]
        h5open(files[i], "r") do f
            # parse out the attributes
            attr = HDF5.attributes(f)
            species[i] = read(attr["species"])
            mass[i] = read(attr["mass"])
            depth[i] = read(attr["depth"])
            λrest[i] = read(attr["air_wavelength"])
            geff[i] = read(attr["g_eff"])
            height[i] = read(attr["height"])
            lower_level[i] = read(attr["lower_level"])
            upper_level[i] = read(attr["upper_level"])
        end
    end
    return LineProperties(species, mass, depth, λrest, geff, height, lower_level, upper_level, filename)
end

get_species(lp::LineProperties) = lp.species
get_rest_wavelength(lp::LineProperties) = lp.λrest
get_energy(lp::LineProperties) = lp.upper_level - lp.lower_level
get_lower_level(lp::LineProperties) = lp.lower_level
get_depth(lp::LineProperties) = lp.depth
get_geff(lp::LineProperties) = lp.geff
get_height(lp::LineProperties) = lp.height
get_file(lp::LineProperties) = lp.file
