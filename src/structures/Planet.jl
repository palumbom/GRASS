struct Planet{T<:AF}
    mass::T
    radius::T
    period::T
    semi_major_axis::T
    eccentricity::T
    inclincation::T
    longitude_ascending_node::T
    longitude_periapsis::T
    mean_longitude::T
end

body_mass(B::Body) = B.mass
T(B::Body) = B.orbital_period
a(B::Body) = B.semi_major_axis
b(B::Body) = B.semi_major_axis * sqrt(1.0 - B.eccentricity^2)
e(B::Body) = B.eccentricity
i(B::Body) = B.inclination
Ω(B::Body) = B.longitude_ascending_node
ωbar(B::Body) = B.longitude_periapsis
L(B::Body) = B.mean_longitude

"""
units:
    - radius: fractional solar radius (rplanet/rstar)
    - period: years
    - semi_major_axis: AU
"""
function Planet(;mass=NaN, radius=NaN, period=NaN, semi_major_axis=NaN, eccentricity=NaN,
                 longitude_ascending_node=NaN, longitude_periapsis=NaN, mean_longitude=NaN)
    @assert !isnan(radius)
    @assert !isnan(period)
    @assert !isnan(semiaxis)
    @assert 0.0 <= inclination <= 180.0

    return Planet(mass, radius, period, semi_major_axis, eccentricity,
                  deg2rad(inclination), deg2rad(longitude_ascending_node),
                  deg2rad(longitude_periapsis), deg2rad(mean_longitude))
end
