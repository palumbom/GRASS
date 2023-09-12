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
    R_y::AA{T,2}
    R_z::AA{T,2}
    O⃗::AA{T,1}
    A::T
    B::T
    C::T
    v0::T
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
function DiskParams(;N=197, Nt=NaN, radius=1.0, inclination=90.0, u1=0.4,
                     u2=0.26, vsini=2067.03346686649251345, A=14.713,
                     B=-2.396, C=-1.787, Nsubgrid=40)
    # assertions and warnings
    @assert !isnan(Nt)

    if N != 197
        @warn "N should be set to 197 for physical validity"
    end

    if Nsubgrid <= 1
        @warn "Nsubgrid must be greater than 1"
        @assert Nsubgrid > 1
    end

    if inclination == 0.0
        @warn "Unresolved bug for i=0, setting i=1e-10"
        inclination = 1e-10
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

    # create rotation matrix for inclination
    @assert -90.0 <= inclination <= 90.0
    iₛ = deg2rad(90.0 - inclination)
    R_x = [1.0 0.0 0.0;
           0.0 cos(iₛ) -sin(iₛ);
           0.0 sin(iₛ) cos(iₛ)]

    # create rotation matrix for stellar rotation
    # TODO consider differential rotation?
    # TODO check math
    ω_eq = vsini / (4.378993e9)
    per = 2π / ω_eq
    v_ang = 2π / per
    dθ = -2π / per
    R_y = [cos(dθ) 0.0 sin(dθ);
           0.0 1.0 0.0;
           -sin(dθ) 0.0 cos(dθ)]


    # create rotation matrix for sky-plane rotation (i.e., position angle)
    pa = deg2rad(0.0)
    R_z = [cos(pa) -sin(pa) 0.0;
           sin(pa) cos(pa) 0.0;
           0.0 0.0 1.0]

    # convert vsini to units of R*/day/c_ms
    v0 = (vsini / c_ms) * (360.0 / A)

    # set observer vector to large distance (units = stellar radius)
    O⃗ = [0.0, 0.0, 1e6]

    return DiskParams(N, Nt, radius, ϕe, ϕc, θe, θc, Nθ, Nsubgrid,
                      R_x, R_y, R_z, O⃗, A, B, C, v0, u1, u2)
end
