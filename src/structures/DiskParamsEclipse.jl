struct DiskParamsEclipse{T<:AF}
    N::Int
    Nt::Int
    ρs::T
    ϕe::AA{T,1}
    ϕc::AA{T,1}
    θe::AA{T,2}
    θc::AA{T,2}
    Nθ::AA{Int,1}
    Nsubgrid::Int
    A::T
    B::T
    C::T
end

"""
    DiskParamsEclipse(; N=197, Nt, Nsubgrid=40, radius=sun_radius,
                      A=14.713, B=-2.396, C=-1.787)

Construct a `DiskParamsEclipse` instance describing the solar disk grid used for
eclipse spectral synthesis. The grid is built in heliographic latitude/longitude;
the observer geometry (and therefore which cells are occulted by the Moon) is
computed per epoch from SPICE ephemerides at synthesis time, rather than being
fixed at construction. This is the eclipse counterpart to [`DiskParams`](@ref).

# Keyword Arguments
- `N=197`: number of stellar latitude grid elements. Should be set to 197 for physical validity; other values emit a warning.
- `Nt`: number of time steps (epochs) to synthesize. **Required** — defaults to a `NaN` sentinel that trips an assertion if not supplied. Must equal `length(time_stamps)` passed to [`synthesize_spectra_eclipse`](@ref).
- `Nsubgrid=40`: subgrid resolution used to refine partially-occulted cells near the lunar limb. Must be greater than 1.
- `radius=sun_radius`: stellar radius in km (defaults to the SPICE solar radius).
- `A=14.713`: differential rotation coefficient (deg/day).
- `B=-2.396`: differential rotation coefficient (deg/day).
- `C=-1.787`: differential rotation coefficient (deg/day).
"""
function DiskParamsEclipse(;N=197, Nt=NaN, Nsubgrid=40, radius=sun_radius,
                            A=14.713, B=-2.396, C=-1.787)
    # assertions and warnings
    @assert !isnan(Nt)

    if N != 197
        @warn "N should be set to 197 for physical validity"
    end

    if Nsubgrid <= 1
        @warn "Nsubgrid must be greater than 1"
        @assert Nsubgrid > 1
    end

    # get latitude grid edges and centers
    ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=N+1)
    ϕc = GRASS.get_grid_centers(ϕe)

    # number of longitudes in each latitude slice
    Nθ = get_Nθ.(ϕc, step(ϕe)) 

    # make longitude grid
    θe = zeros(N+1, maximum(Nθ)+1)
    θc = zeros(N, maximum(Nθ))
    for i in eachindex(Nθ)
        edges = range(deg2rad(0.0), deg2rad(360.0), length=Nθ[i]+1)
        θc[i, 1:Nθ[i]] .= GRASS.get_grid_centers(edges)
        θe[i, 1:Nθ[i]+1] .= collect(edges)
    end

    return DiskParamsEclipse(N, Nt, radius, ϕe, ϕc, θe, θc, Nθ, Nsubgrid, A, B, C)
end

function get_Nθ(ϕc, dϕ)
    return ceil(Int, 2π * cos(ϕc) / dϕ)
end