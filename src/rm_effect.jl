struct Planet{T<:AF}
    radius::T
    period::T
    semiaxis::T
    eccentricity::T
    inclincation::T
    vcirc::T
    b::T
end

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

function is_transiting(p::Planet{T}) where T<:AF
    return abs(p.b) <= 1.0
end

function plot_planet_transit(disk::DiskParams, planet::Planet)
    # get grids
    grid_1D = make_grid(disk.N)
    grid_2D = make_grid_2D(grid_1D)
    grid_xs = get_grid_xs(grid_2D)
    grid_ys = get_grid_ys(grid_2D)
    grid_edges = get_grid_edges(grid_1D)

    # prodi = Iterators.product(grid, grid)
    # grid_xs = get_grid_xs(prodi)
    # grid_ys = get_grid_ys(prodi)


    # # get map of rotational velocities
    # mus = calc_mu.(prodi)
    # vels = patch_velocity_los.(prodi, pole=disk.pole) .* c_kms
    # plt.imshow(vels); plt.show()
    # ints = calc_norm_terms(disk)
    # ints[mus .<= 0.0] .= NaN
    # vels[mus .<= 0.0] .= NaN

    # println("derp")

    # # customize color map
    # cmap = pycopy.copy(mpl.pyplot.matplotlib.cm.get_cmap("seismic"))
    # cmap.set_bad("k")
    # # cmap.set_under("k")

    # # initialize fig object
    # fig = plt.figure()
    # axs = fig.add_subplot()

    # # plot the disk with line for planet
    # # img = axs.imshow(vels', origin="lower", cmap=cmap,
    # #                  extent=[-1,1,-1,1],
    # #                  vmin=minimum(filter(!isnan, vels)),
    # #                  vmax=maximum(filter(!isnan, vels)))
    # img = axs.pcolormesh(grid_xs, grid_ys, vels, cmap=cmap, shading="auto",
    #                      vmin=minimum(filter(!isnan, vels)),
    #                      vmax=maximum(filter(!isnan, vels)))
    # axs.axhline(planet.b, c="k")
    # cb = fig.colorbar(img)
    # cb.ax.set_ylabel("LOS Velocity (km/s)")
    # plt.show()
    return nothing
end
