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
- `N=197`: Number of stellar latitude grid elements. Should be set to 197 for physical validity.
- `Nt=50`: Number of 15-second time steps.
- `radius=1.0`: Radius of model star. Default is one solar radius.
- `inclination=90.0`: Sky-plane inclination of model stellar disk. 90 degrees is equator on.
- `u1=0.4`: Quadratic limb darkening law coefficient.
- `u2=0.25`: Quadratic limb darkening law coefficient
- `vsini=2067.03346686649251345`: Equatorial rotational velocity magnitude in units of meters per second.
- `A=14.713`: Differential rotation coefficient. Units of deg/day.
- `B=-2.396`: Differential rotation coefficient. Units of deg/day.
- `C=-1.787`: Differential rotation coefficient. Units of deg/day.
- `dist=4.435e7`: Distance to observer. Default is one parsec in solar radii.
"""
function DiskParams(;N=197, Nt=NaN, Nsubgrid=40, radius=sun_radius,
                     u1=0.4, u2=0.26, A=14.713, B=-2.396, C=-1.787)
    # assertions and warnings
    @assert !isnan(Nt)

    if N != 197
        @warn "N should be set to 197 for physical validity"
    end

    if Nsubgrid <= 1
        @warn "Nsubgrid must be greater than 1"
        @assert Nsubgrid > 1
    end

    # if inclination == 0.0
    #     @warn "Unresolved bug for i=0, setting i=1e-10"
    #     inclination = 1e-10
    # end

    # get latitude grid edges and centers
    ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=N+1)
    ϕc = get_grid_centers(ϕe)

    # number of longitudes in each latitude slice
    Nθ = get_Nθ.(ϕc, step(ϕe)) #ceil.(Int, 2π .* cos.(ϕc) ./ step(ϕe))

    # make longitude grid
    θe = zeros(N+1, maximum(Nθ)+1)
    θc = zeros(N, maximum(Nθ))
    for i in eachindex(Nθ)
        edges = range(deg2rad(0.0), deg2rad(360.0), length=Nθ[i]+1)
        θc[i, 1:Nθ[i]] .= get_grid_centers(edges)
        θe[i, 1:Nθ[i]+1] .= collect(edges)
    end

    # # create rotation matrix for inclination
    # @assert -90.0 <= inclination <= 90.0
    # iₛ = deg2rad(90.0 - inclination)
    # R_x = [1.0 0.0 0.0;
    #        0.0 cos(iₛ) -sin(iₛ);
    #        0.0 sin(iₛ) cos(iₛ)]

    # # create rotation matrix for stellar rotation
    # # TODO consider differential rotation?
    # # TODO check math
    # ω_eq = vsini / (4.378993e9)
    # per = 2π / ω_eq
    # v_ang = 2π / per
    # dθ = -2π / per
    # R_y = [cos(dθ) 0.0 sin(dθ);
    #        0.0 1.0 0.0;
    #        -sin(dθ) 0.0 cos(dθ)]


    # # create rotation matrix for sky-plane rotation (i.e., position angle)
    # pa = deg2rad(0.0)
    # R_z = [cos(pa) -sin(pa) 0.0;
    #        sin(pa) cos(pa) 0.0;
    #        0.0 0.0 1.0]

    # # convert vsini to units of R*/day/c_ms
    # v0 = (vsini / c_ms) * (360.0 / A)

    # # set vector pointing from observer to star center
    # O⃗ = [0.0, 0.0, dist]

    # return DiskParams(N, Nt, radius, ϕe, ϕc, θe, θc, Nθ, Nsubgrid,
    #                   R_x, R_y, R_z, O⃗, A, B, C, v0, u1, u2)
    return DiskParams(N, Nt, radius, ϕe, ϕc, θe, θc, Nθ, Nsubgrid, A, B, C, u1, u2)
end

function get_Nθ(ϕc, dϕ)
    return ceil(Int, 2π * cos(ϕc) / dϕ)
end
