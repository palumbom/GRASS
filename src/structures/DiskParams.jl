struct DiskParams{T<:AF}
    N::Int
    Nt::Int
    ρs::T
    ϕe::AA{T,1}
    ϕc::AA{T,1}
    θe::AA{T,2}
    θc::AA{T,2}
    Nθ::AA{Int,1}
    Nsubgrid::Int
    R_x::AA{T,2}
    O⃗::AA{T,1}
    A::T
    B::T
    C::T
    u1::T
    u2::T
end

"""
    DiskParams(; N=132, Nt=50, inclination=90.0)

Construct a `DiskParams` composite type instance. In the coordinate system,
the x- and z- axes are sky-plane, and the y-axis is along the observer-to-star-
center vector.

# Arguments
- `N::Integer=132`: Number of stellar latitude grid elements
- `Nt::Integer=50`: Number of 15-second snapshots.
- `Inclination::Float64`: Sky-plane inclination of stellar grid. 90.0 is equator-on.
"""
function DiskParams(;N=275, Nt=NaN, radius=1.0, inclination=90.0, u1=0.4,
                     u2=0.26, A=14.713, B=-2.396, C=-1.787, Nsubgrid=40)
    # assertions and warnings
    @assert !isnan(Nt)
    @assert Nsubgrid > 1

    if N != 275
        @warn "N should be set to 275 for physical validity!"
    end

    # get latitude grid edges and centers
    ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=N+1)
    ϕc = get_grid_centers(ϕe)

    # number of longitudes in each latitude slice
    Nθ = ceil.(Int, 2π .* cos.(ϕc) ./ step(ϕe))

    # make longitude grid
    θe = zeros(N+1, maximum(Nθ)+1)
    θc = zeros(N, maximum(Nθ))
    for i in eachindex(Nθ)
        edges = range(deg2rad(0.0), deg2rad(360.0), length=Nθ[i]+1)
        θc[i, 1:Nθ[i]] .= get_grid_centers(edges)
        θe[i, 1:Nθ[i]+1] .= collect(edges)
    end

    # create rotation matrix
    @assert -90.0 <= inclination <= 90.0
    iₛ = deg2rad(90.0 - inclination)
    R_x = M = [1.0 0.0 0.0;
               0.0 cos(iₛ) sin(iₛ);
               0.0 -sin(iₛ) cos(iₛ)]

    # set observer vector to large distance (units = stellar radius)
    # O⃗ = [0.0, 220.0, 0.0] # 220 => 1 AU in Solar Radii
    O⃗ = [0.0, 1e6, 0.0]

    return DiskParams(N, Nt, radius, ϕe, ϕc, θe, θc, Nθ, Nsubgrid, R_x, O⃗, A, B, C, u1, u2)
end
