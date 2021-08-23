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

function quad_limb_darkening(x::T, y::T, u1::T, u2::T) where T<:AF
    return quad_limb_darkening(calc_mu(x,y), u1, u2)
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
function rotation_period(sin_lat::T; A::T=14.713, B::T=-2.396, C::T=-1.787) where T<:AF
    @assert -1.0 <= sin_lat <= 1.0
    return 360.0/(A + B * sin_lat^2 + C * sin_lat^4)
end

"""
    patch_velocity_los(x, y; rstar, pole)

Compute line of sight velocity of a patch of stellar surface given by x,y (assumed in [-1,1]).
Return value is in (Rsol/day)/speed_of_light (i.e., dimensionless like z = v/c)

# Arguments
- `rstar::Float64=1.0`: in R_sol (affects return velocity, but not x,y)
- `pole = (0,1,0)`: unit vector for stellar rotation axis (default is equator-on)
"""
function patch_velocity_los(x::T,y::T; rstar=one(T), pole=(zero(T), one(T), zero(T))) where T<:AF
    polex, poley, polez = pole
    v0 = 0.000168710673 # in (Rsol/day)/speed_of_light
    z = sqrt(rstar - calc_r2(x,y))
    sin_lat = (x*polex) + (y*poley) + (z*polez)
    vmax = v0*rstar/rotation_period(sin_lat)
    return vmax * (polex*y - poley*x)
end

function patch_velocity_los(t::Tuple{T,T}; rstar=one(T), pole=(zero(T), one(T), zero(T))) where T<:AF
    return patch_velocity_los(t[1], t[2], rstar=rstar, pole=pole)
end
