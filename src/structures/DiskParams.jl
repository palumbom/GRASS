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
    R_θ::AA{T,2}
    O⃗::AA{T,1}
    A::T
    B::T
    C::T
    u1::T
    u2::T
end

"""
    DiskParams(; N=132, Nt=50, pole=(0.0, 1.0, 0.0))

Construct a `DiskParams` composite type instance.

# Arguments
- `N::Integer=132`: Length of N*N spatial grid
- `Nt::Integer=50`: Number of 15-second snapshots
- `pole::Tuple{Float64, Float64, Float64}=(0.0, 1.0, 0.0)`: Unit vector specificying rotation axis direction. Default is equator-on.
"""
function DiskParams(;N=132, Nt=NaN, radius=1.0, inclination=90.0, u1=0.4,
                     u2=0.26, A=14.713, B=-2.396, C=-1.787, Nsubgrid=50)
    # assertions and warnings
    @assert !isnan(Nt)
    # if N != 132
    #     @warn "N should be set to 132 for physical validity!"
    # end

    # get latitude grid edges and centers
    ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=N+1)
    ϕc = get_grid_centers(ϕe)

    # calculate span of tile in latitude in fraction of stellar radius
    w_tile = 2.0 * π * radius / N

    # number of longitudes in each latitude slice
    Nθ = ceil.(Int, 2.0 .* π .* radius .* cos.(ϕc) ./ w_tile)

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
    R_θ = M = [1.0 0.0 0.0;
               0.0 cos(iₛ) sin(iₛ);
               0.0 -sin(iₛ) cos(iₛ)]

    # set observer vector to large distance (units = stellar radius)
    # O⃗ = [0.0, 220.0, 0.0] # 220 => 1 AU in Solar Radii
    O⃗ = [0.0, 1e6, 0.0]

    return DiskParams(N, Nt, radius, ϕe, ϕc, θe, θc, Nθ, Nsubgrid, R_θ, O⃗, A, B, C, u1, u2)
end
