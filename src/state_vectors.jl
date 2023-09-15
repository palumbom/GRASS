"""
Calculates cartesian state vector from orbital elements and eccentric anomaly.
See https://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node34.html

# Arguments
- `body::Body`: an instance of the body custom type
- `epoch::Float64`: time elapsed since periastron reference epoch

# Outputs
- `XYZ::Array{Float64,1}`: 3-vector of initial cartesian coordinates
- `XYZ_dot::Array{Float64,1}`: 3-vector of initial cartesian velocities
"""
function calc_state_vector!(xyz::AA{T1,1}, xyz_dot::AA{T1,1}, body::Planet{T1}, epoch::T1) where T1<:AF
    # get the reduced mass of the body
    μ = mass(body) / (mass(body) + 1.0)

    # store semi-major and -minor axes in memory
    major = a(body)

    # get mean motion
    if isnan(body.period)
        n = 2π * sqrt(μ / major^3.0)
    else
        n = 2π / period(body)
    end

    # get eccentricity
    ecc = e(body)

    # get longitude of periapsis
    ω̄ = Ω(body) + ω(body)

    # get the mean anomaly at epoch
    M0 = L(body) - ω̄
    M = mod(M0 + n * epoch, 2π)

    # solve kepler's equation for eccentric anomaly
    E = calc_ecc_anom_iterative_laguerre(M, ecc)

    # enforce that eccentric anomaly converged
    @assert isapprox(kepler(E, ecc, M), 0.0, atol=1e-10)

    # get the true anomaly
    ν = 2.0 * atan(sqrt(1.0 + ecc) * sin(E / 2.0),
                   sqrt(1.0 - ecc) * cos(E / 2.0))

    # get distance to central body
    rt = major * (1.0 - ecc * cos(E))

    # convert to solar radii
    rt *= 1.496e11/ 6.957e8

    # assign position in body frame
    xyz[1] = rt .* cos(ν)
    xyz[2] = rt .* sin(ν)
    xyz[3] = 0.0

    # get velocity magnitude
    # TODO I have no idea if this is right or what units these are
    vt = sqrt(μ * major) / rt

    # assign velocity vector
    xyz_dot[1] = vt * (-sin(E))
    xyz_dot[2] = vt * (sqrt(1.0 - ecc^2.0) * cos(E))
    xyz_dot[3] = 0.0

    # shortcut trig evals
    cosω = cos(ω(body))
    sinω = sin(ω(body))
    cosi = cos(i(body))
    sini = sin(i(body))
    cosΩ = cos(Ω(body))
    sinΩ = sin(Ω(body))

    # copy values before rotation
    x0 = deepcopy(xyz[1])
    y0 = deepcopy(xyz[2])
    z0 = deepcopy(xyz[3])

    # rotate into observer frame
    xyz[1] = x0 * (cosω * cosΩ - sinω * cosi * sinΩ) - y0 * (sinω * cosΩ + cosω * cosi * sinΩ)
    xyz[2] = x0 * (cosω * sinΩ + sinω * cosi * cosΩ) + y0 * (cosω * cosi * cosΩ - sinω * sinΩ)
    xyz[3] = x0 * (sinω * sini) + y0 * (cosω * sini)

    # copy values before rotation
    x0 = deepcopy(xyz_dot[1])
    y0 = deepcopy(xyz_dot[2])
    z0 = deepcopy(xyz_dot[3])

    # rotate into observer frame
    xyz_dot[1] = x0 * (cosω * cosΩ - sinω * cosi * sinΩ) - y0 * (sinω * cosΩ + cosω * cosi * sinΩ)
    xyz_dot[2] = x0 * (cosω * sinΩ + sinω * cosi * cosΩ) + y0 * (cosω * cosi * cosΩ - sinω * sinΩ)
    xyz_dot[3] = x0 * (sinω * sini) + y0 * (cosω * sini)

    return nothing
end

function calc_state_vector!(ros_allocs::RossiterAllocs{T1}, body::Planet{T1}) where T1<:AF
    # alias rossiter allocs
    xyz_planet = ros_allocs.xyz_planet
    xyz_dot_planet = ros_allocs.xyz_dot_planet
    xyz_star = ros_allocs.xyz_star
    xyz_dot_star = ros_allocs.xyz_dot_star
    epochs = ros_allocs.epochs

    # get epochs
    epochs .= collect(range(0.0, length=length(epochs), step = 15.0 / 3.154e7))

    # loop over epoch
    for t in eachindex(epochs)
        calc_state_vector!(view(xyz_planet, :, t), view(xyz_dot_planet, :, t), body, epochs[t])
    end
    return nothing
end
