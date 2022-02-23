struct Planet{T<:AF}
    radius::T
    period::T
    semiaxis::T
    eccentricity::T
    inclincation::T
    vcirc::T
    b::T
end

const rsun_au = 0.00465
const sec_year = 3.154e7

"""
units:
    - radius: fractional solar radius (rplanet/rstar)
    - period: years
    - semiaxis: AU
"""
function Planet(;radius=NaN, period=NaN, semiaxis=NaN, inclination=90.0)
    @assert !isnan(radius)
    @assert !isnan(period)
    @assert !isnan(semiaxis)
    @assert 0.0 <= inclination <= 180.0

    # get circular velocity from period, convert to solar radii/s
    vcirc = 2π * semiaxis / period
    vcirc *= 1.0 / (rsun_au * sec_year)

    # calculate impact parameter
    b = (semiaxis / rsun_au) * cos(deg2rad(inclination))
    if abs(b) > 1.0
        println(">>> Planet will not transit!")
    end
    return Planet(radius, period, semiaxis, 0.0, inclination, vcirc, b)
end

function calc_planet_position(t::T, planet::Planet{T}) where T<:AF
    return planet.vcirc * (t) - 1.0, planet.b
end

function calc_planet_position(t::AA{T,1}, planet::Planet{T}) where T<:AF
    out = map(x -> calc_planet_position(x, planet), t)
    return map(collect, zip(out...))
end

function calc_transit_duration(p::Planet{T}) where T<:AF
    t1 = p.period / (π * p.semiaxis)
    t2 = sqrt((1.0 + p.radius)^2 - (p.semiaxis * cos(p.inclination)^2))
    return t1 * t2
end
