struct DiskParams{T<:AF}
    N::Int
    Nt::Int
    ρs::T
    ϕe::AA{T,1}
    θe::AA{T,1}
    ϕc::AA{T,1}
    θc::AA{T,1}
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
function DiskParams(;N=132, Nt=NaN, inclination=90.0, u1=0.4,
                     u2=0.26, A=14.713, B=-2.396, C=-1.787)
    # assertions and warnings
    @assert !isnan(Nt)
    # if N != 132
    #     @warn "N should be set to 132 for physical validity!"
    # end

    # get grid edges
    ϕe, θe = make_grid(N)

    # get grid centers
    ϕc = get_grid_centers(ϕe)
    θc = get_grid_centers(θe)

    # create rotation matrix
    @assert -90.0 <= inclination <= 90.0
    iₛ = deg2rad(90.0 - inclination)
    R_θ = M = [1.0 0.0 0.0;
               0.0 cos(iₛ) sin(iₛ);
               0.0 -sin(iₛ) cos(iₛ)]

    # set observer vector to Earth-Sun distance in AU
    O⃗ = [0.0, 220.0, 0.0]

    return DiskParams(N, Nt, 1.0, ϕe, θe, ϕc, θc, R_θ, O⃗, A, B, C, u1, u2)
end
