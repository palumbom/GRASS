"""
Author: Michael Palumbo
Created: April 2019
Contact: mlp95@psu.edu
"""

"""
Absorption line
"""
function absorption_line(x::T; mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    return one(T) - depth * exp(-((x-mid)/width)^2.0/2.0)
end

function lorentzian_line(x::T; mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    y = (x - mid)/(width/2)
    return one(T) - (depth / (one(T) + y^2))
end

"""

"""
function doppler_factor(vel::T) where T<:AF
    num = one(T) + vel/c_ms
    den = one(T) - vel/c_ms
    return sqrt(num/den)
end

"""
Width of line due to thermal broadening
"""
function width_thermal(; λ::T1=1.0, M::T1=1.0, T::T1=5778.0, v_turb::T1=0.0) where T1<:AF
    return sqrt((2.0*kB/mH) * (T/M) + v_turb^2) * (λ/c)
end

"""
Quadratic limb darkening law.
Takes μ = cos(heliocentric angle) and LD parameters, u1 and u2.
"""
function quad_limb_darkening(μ::T, u1::T, u2::T) where T<:AF
    return !iszero(μ) * (one(T) - u1*(one(T)-μ) - u2*(one(T)-μ)^2)
end

function quad_limb_darkening(x::T, y::T, u1::T, u2::T) where T<:AF
    return quad_limb_darkening(calc_mu(x,y), u1, u2)
end

"""
Integrate over smooth stellar disk assuming quadratic limb darkening law.
Uses N^2 evaluations over square.
"""
function calc_intensity_normalization(; u1::T=0.4, u2::T=0.26, N::Int64=256) where T<:AF
    # make grid on which to calculate intensities
    grid = range(-1, 1, length=N)
    f(t::Tuple) = quad_limb_darkening(t[1], t[2], u1, u2)

    # find sum of intensities, normalization factor
    sumI = mapreduce(f, +, ((x,y) for x in grid, y in grid))
    sumI *= pi / length(grid)^2
    return sumI
end

"""
Calculate stellar rotation period in days:
- sin_lat = sine of stellar latitude
- A, B, C = Diff. rotation equation parameters (in deg/day)
"""
function rotation_period(sin_lat::T; A::T=14.713, B::T=-2.396, C::T=-1.787) where T<:AF
    @assert -1.0 <= sin_lat <= 1.0
    return 360.0/(A + B * sin_lat^2 + C * sin_lat^4)
end

"""
patch_velocity_los(x, y; rstar, pole)
Compute line of sight velocity of a patch of stellar surface given by x,y (assumed in [-1,1]).
Optional params:
- rstar = 1: in R_sol (affects return velocity, but not x,y)
- pole = (0,1,0): unit vector with stellar rotation axis
Depends on:
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

"""

"""
function wav2vel(λ::T, λ0::T) where T<:AF
    return (λ - λ0) * c_ms / λ0
end
