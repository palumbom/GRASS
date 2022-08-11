abstract type AbstractLineProperties end
struct LineProperties{T<:AF} <: AbstractLineProperties
    species::String
    mass::T             # amu
    depth::T
    λrest::T            # angstroms
    geff::T
    height::T           # km
    lower_level::T      # eV
    upper_level::T      # eV
end

function LineProperties(filename::String)
    @assert isfile(filename)
    lp = h5open(filename, "r") do f
        g = f["properties"]
        attr = HDF5.attributes(g)
        species = read(attr["species"])
        mass = read(attr["mass"])
        depth = read(attr["depth"])
        λrest = read(attr["air_wav"])
        geff = read(attr["g_eff"])
        height = read(attr["height"])
        lower_level = read(attr["lower_level"])
        upper_level = read(attr["upper_level"])
        return LineProperties(species, mass, depth, λrest, geff, height, lower_level, upper_level)
    end
    return lp
end

get_species(lp::LineProperties) = lp.species
get_rest_wavelength(lp::LineProperties) = lp.λrest
get_energy(lp::LineProperties) = lp.upper_level - lp.lower_level
get_lower_level(lp::LineProperties) = lp.lower_level
get_depth(lp::LineProperties) = lp.depth
get_geff(lp::LineProperties) = lp.geff
get_height(lp::LineProperties) = lp.height
