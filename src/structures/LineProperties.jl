abstract type AbstractLineProperties end
struct LineProperties{T<:AF} <: AbstractLineProperties
    species::String
    mass::T             # amu
    depth::T
    λrest::T            # angstroms
    geff::T
    height::T           # km
end
