abstract type AbstractLineProperties end
struct LineProperties{T<:AF} <: AbstractLineProperties
    species::String
    mass::T             # amu
    depth::T
    λrest::T            # angstroms
    geff::T
    height::T           # km
end

function LineProperties(filename::String)
    @assert isfile(filename)

    lp = h5open(filename, "r") do f
        g = f["properties"]
        attr = HDF5.attributes(g)
        species = read(attr["species"])
        mass = read(attr["mass"])
        depth = read(attr["depth"])
        λrest = read(attr["air_wavelength"])
        geff = read(attr["g_eff"])
        height = read(attr["height"])
        return LineProperties(species, mass, depth, λrest, geff, height)
    end
    return lp
end

get_rest_wavelength(lp::LineProperties) = lp.λrest
