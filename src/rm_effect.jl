struct Planet{T<:AF}
    radius::T
    period::T
    semiaxis::T
    vcirc::T
    eccentricity::T
    pos::MVector{2,T}
end

"""
units:
    - radius: fractional solar radius (rplanet/rstar)
    - period: years
    - semiaxis: AU
"""
function Planet(;radius=NaN, period=NaN, semiaxis=NaN)
    # get circular velocity from period, convert to solar radii/s
    vcirc = 2π * semiaxis / period
    vcirc *= (0.00465/3.154e7)

    # set initial position
    pos = @MVector [-1.0, 0.0]
    return Planet(radius, period, semiaxis, vcirc, 0.0, pos)
end

function calc_planet_position(t::T, planet::Planet{T}) where T<:AF

    return x_planet, y_planet
end

function is_occulted(x_star::T, y_star::T, x_planet::T,
                     y_planet::T, r_planet::T) where T<:AF
    r2 = calc_r2(x_star, y_star)
    dist = sqrt((x_star - x_planet)^2 + (y_star - y_planet)^2)
    return (r2 < 1.0) * (dist < r_planet)
end
