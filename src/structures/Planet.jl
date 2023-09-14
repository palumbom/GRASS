struct Planet{T1<:AF}
    mass::T1
    radius::T1
    period::T1
    semi_major_axis::T1
    eccentricity::T1
    inclination::T1
    longitude_ascending_node::T1
    argument_periapsis::T1
    mean_longitude::T1
end

mass(p::Planet) = p.mass
period(p::Planet) = p.period
a(p::Planet) = p.semi_major_axis
b(p::Planet) = p.semi_major_axis * sqrt(1.0 - p.eccentricity^2)
e(p::Planet) = p.eccentricity
i(p::Planet) = p.inclination
Ω(p::Planet) = p.longitude_ascending_node
ω(p::Planet) = p.argument_periapsis
L(p::Planet) = p.mean_longitude

"""
units:
    - mass: mass in solar masses
    - radius: planet radius as fraction of stellar radius
    - period: period in years
    - semi_major_axis: semi major axis in AU
    - eccentricty
    - inclination
    - longitude_ascending_node
    - longitude_periapsis
    - mean_longitude
"""
function Planet(;mass=NaN, radius=NaN, period=NaN, semi_major_axis=NaN,
                 eccentricity=NaN, inclination=NaN, longitude_ascending_node=NaN,
                 longitude_periapsis=NaN, mean_longitude=NaN)
    @assert !isnan(radius)
    @assert !isnan(period)
    @assert !isnan(semi_major_axis)
    @assert 0.0 <= inclination <= 180.0

    return Planet(mass, radius, period, semi_major_axis, eccentricity,
                  deg2rad(inclination), deg2rad(longitude_ascending_node),
                  deg2rad(longitude_periapsis), deg2rad(mean_longitude))
end
