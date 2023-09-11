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
function calc_state_vector(body::Planet{T1}; epoch=0.0) where T1<:AF
    # calculate argument of periapsis and store eccentricity in memory
    ω = ω̄(body) - Ω(body)
    ecc = e(body)

    # get mean motion
    n = 2π / T(body)

    # get mean anomaly
    M = mod(n * (epoch), 2π) # L_t - ω̄(body)

    # iteratively solve Kepler's equation for eccentric anomaly
    E = calc_ecc_anom_iterative_laguerre(M, ecc)

    # enforce that eccentric anomaly converged
    @assert isapprox(kepler(E, ecc, M), 0.0, atol=1e-10)

    # store trig evals in memory
    cosE = cos(E)
    sinE = sin(E)
    cosω = cos(ω)
    sinω = sin(ω)
    cosi = cos(i(body))
    sini = sin(i(body))
    cosΩ = cos(Ω(body))
    sinΩ = sin(Ω(body))

    # store semi-major and -minor axes in memory
    major = a(body)
    minor = b(body)

    # calclate position in 2D orbital plane
    xyz = [major * (cosE - ecc), minor * sinE, 0]

    # construct rotation matrices
    Ω_mat = [cosΩ -sinΩ 0.0; sinΩ cosΩ 0.0; 0.0 0.0 1.0]
    i_mat = [1.0 0.0 0.0; 0.0 cosi -sini; 0.0 sini cosi]
    ω_mat = [cosω -sinω 0.0; sinω cosω 0.0; 0.0 0.0 1.0]

    rot_mat = Ω_mat * i_mat * ω_mat

    # perform rotation into 3D vector
    XYZ = rot_mat * xyz

    if a(body) == 0
        return XYZ, [0.0, 0.0, 0.0]
    else
        # rate of change of eccentric anomaly
        Mdot = 1.0 / T(body)
        Edot = Mdot / (1.0 - ecc * cosE)

        # calculate velocity in 2D orbital plane
        xyz_dot = [-major * sinE * Edot, minor * cosE * Edot, 0]

        # perform rotation into 3D vector
        XYZ_dot = rot_mat * xyz_dot
    end

    return XYZ, XYZ_dot
end


# function calc_planet_position(t::T, planet::Planet{T}) where T<:AF
#     return planet.vcirc * (t) - 1.0, planet.b
# end

# function calc_planet_position(t::AA{T,1}, planet::Planet{T}) where T<:AF
#     out = map(x -> calc_planet_position(x, planet), t)
#     return map(collect, zip(out...))
# end

# function calc_transit_duration(p::Planet{T}) where T<:AF
#     t1 = p.period / (π * p.semiaxis)
#     t2 = sqrt((1.0 + p.radius)^2 - (p.semiaxis * cos(p.inclination)^2))
#     return t1 * t2
# end

# function is_transiting(p::Planet{T}) where T<:AF
#     return abs(p.b) <= 1.0
# end

# function plot_planet_transit(disk::DiskParams, planet::Planet)
#     # get grids
#     grid_1D = make_grid(disk.N)
#     grid_2D = make_grid_2D(disk.N)
#     grid_xs = get_grid_xs(grid_2D)
#     grid_ys = get_grid_ys(grid_2D)
#     grid_edges = get_grid_edges(grid_1D)

#     # get map of rotational velocities
#     mus = calc_mu.(grid_2D)
#     vels = patch_velocity_los.(grid_2D, pole=disk.pole) .* c_kms
#     ints = calc_norm_terms(disk)

#     ints[mus .<= 0.0] .= NaN
#     vels[mus .<= 0.0] .= NaN

#     # customize color map
#     cmap = pycopy.copy(mpl.pyplot.matplotlib.cm.get_cmap("seismic"))
#     cmap.set_bad("k")
#     # cmap.set_under("k")

#     # initialize fig object
#     fig = plt.figure()
#     axs = fig.add_subplot()

#     # plot the disk with line for planet
#     img = axs.pcolormesh(grid_xs, grid_ys, vels, cmap=cmap, shading="auto",
#                          vmin=minimum(filter(!isnan, vels)),
#                          vmax=maximum(filter(!isnan, vels)))
#     axs.axhline(planet.b, c="k")
#     cb = fig.colorbar(img)
#     cb.ax.set_ylabel("LOS Velocity (km/s)")
#     plt.show()
#     return nothing
# end
