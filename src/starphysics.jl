function absorption_line(x::T; mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    return one(T) - depth * exp(-((x-mid)/width)^2.0/2.0)
end

function gaussian_line(x::T; mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    return absorption_line(x, mid=mid, width=width, depth=depth)
end

function lorentzian_line(x::T; mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    y = (x - mid)/(width/2)
    return one(T) - (depth / (one(T) + y^2))
end

# calculate doppler factor
function doppler_factor(vel::T) where T<:AF
    num = one(T) + vel/c_ms
    den = one(T) - vel/c_ms
    return sqrt(num/den)
end

# Width of line due to thermal broadening
function width_thermal(; λ::T1=1.0, M::T1=1.0, T::T1=5778.0, v_turb::T1=0.0) where T1<:AF
    return sqrt((2.0*kB/mH) * (T/M) + v_turb^2) * (λ/c)
end

# Quadratic limb darkening law.
# Takes μ = cos(heliocentric angle) and LD parameters, u1 and u2.
function quad_limb_darkening(μ::T, u1::T, u2::T) where T<:AF
    return !iszero(μ) * (one(T) - u1*(one(T)-μ) - u2*(one(T)-μ)^2)
end

"""
    rotation_period(sin_lat; A=14.713, B=2.396, C=1.787)

Calculate stellar rotation period in days at given sine of latitude.

# Arguments
- `sin_lat::Float64`: Sine of stellar latitude.
- `A::Float64=14.713`: Coefficient from Snodgrass & Ulrich (1990)
- `B::Float64=2.396`: Coefficient from Snodgrass & Ulrich (1990)
- `C::Float64=1.787`: Coefficient from Snodgrass & Ulrich (1990)
"""
function rotation_period(ϕ::T; A::T=14.713, B::T=-2.396, C::T=-1.787) where T<:AF
    @assert -π/2 <= ϕ <= π/2
    sinϕ = sin(ϕ)
    return 360.0/(A + B * sinϕ^2.0 + C * sinϕ^4.0)
end

"""
    patch_velocity_los(x, y; rstar, pole)

Compute line of sight velocity of a patch of stellar surface given by x,y (assumed in [-1,1]).
Return value is in (Rsol/day)/speed of light (i.e., dimensionless like z = v/c)

# Arguments
- `rstar::Float64=1.0`: in R_sol (affects return velocity, but not x,y)
- `pole = (0,1,0)`: unit vector for stellar rotation axis (default is equator-on)
"""
function patch_velocity_los(ϕ::T, θ::T, disk::DiskParams{T}; P⃗=[0.0, 0.0, disk.ρs]) where T<:AF
    # get vector pointing from spherical circle to patch
    xyz = sphere_to_cart(disk.ρs, ϕ, θ)

    # get vector from spherical circle center to surface patch
    C⃗ = xyz .- [0.0, 0.0, last(xyz)]

    # velocity magnitude at equator, in Rsol/day/c_ms
    v0 = 0.000168710673

    # get velocity vector direction and set magnitude
    vel = cross(C⃗, P⃗)
    vel /= norm(vel)
    vel *= (v0 / rotation_period(ϕ; A=disk.A, B=disk.B, C=disk.C))

    # rotate by stellar inclination
    xyz .= disk.R_θ * xyz
    vel .= disk.R_θ * vel

    # find get vector from observer to surface patch, return projection
    O⃗_surf = xyz .- disk.O⃗
    angle = dot(O⃗_surf, vel) / (norm(O⃗_surf) * norm(vel))
    return norm(vel) * angle
end
