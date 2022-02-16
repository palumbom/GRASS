struct Planet{T<:AF}
    radius::T
    period::T
end

function calc_planet_position(t::T, planet::Planet{T}) where T<:AF


    return x_planet, y_planet
end

function is_occulted(x_star::T, y_star::T, x_planet::T,
                     y_planet::T, r_planet::T) where T<:AF
    r2 = calc_r2(x_star, y_star)
    dist = sqrt((x_star - x_planet)^2 + (y_star - y_planet)^2)
    return (r2 < 1.0) * (dist > r_planet)
end
