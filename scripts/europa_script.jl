# load packages
# using Pkg; Pkg.activate(".")
using Glob
using GRASS
using SPICE
using PyCall
using NaNMath
using Statistics
using LaTeXStrings
using LinearAlgebra
import PyPlot
plt = PyPlot
mpl = plt.matplotlib
import Base.Iterators

rcParams = PyDict(mpl["rcParams"])
rcParams["text.usetex"] = true

# define scattering function
function lambertian(mu0::T) where T<:Real
    return mu0
end

function lommel_seeliger(mu0::T, mu::T) where T<:Real
    return 2 * mu0 / (mu0 + mu)
end

# get the kernels
GRASS.get_kernels()

# pull out constants for body dimensions
earth_radius = bodvrd("EARTH", "RADII")[1]
earth_radius_pole = bodvrd("EARTH", "RADII")[3]
sun_radius = bodvrd("SUN","RADII")[1]
moon_radius = bodvrd("MOON", "RADII")[1]
europa_radius = bodvrd("EUROPA", "RADII")[1]

# set string for ltt and abberation
lt_flag = "CN+S"
# lt_flag = "NONE"

# set HARPS-N location
obs_lat = 28.754
obs_long = -17.88814
alt = 2.370

# get state vector of observatory in Earth frame
flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)
EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])

# set obstime
tstart = "2014-01-05T14:00:00.00"
tstop = "2014-01-06T04:00:00.00"
obstimes = range(utc2et(tstart), utc2et(tstop), length=24)

t_ref = "2014-01-05T23:50:00.00"
obstime_ref = utc2et(t_ref)

# set up paramaters for solar disk
N = 50
Nt = length(obstimes)
Nsubgrid = 4
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=Nsubgrid)

# set up workspace
wsp = GRASS.SynthWorkspaceEclipse(disk, 1, Nt, verbose=true)
ϕc = wsp.ϕc
θc = wsp.θc
μs = wsp.μs
ld = wsp.ld
dA = wsp.dA
wts = wsp.wts
z_rot = wsp.z_rot
mean_weight_v_no_cb = similar(z_rot)
scattering_lambert = similar(z_rot)
scattering_super_lambert = similar(z_rot)

# allocate memory for vectors
SP_sun_pos = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
SP_sun_vel = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
SP_bary = fill(Vector{Float64}(undef, 6), Nsubgrid, Nsubgrid)
SP_bary_pos = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
SP_bary_vel = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
pole_vector_grid = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
J2P_bary = fill(Vector{Float64}(undef, 6), Nsubgrid, Nsubgrid)

# allocate memory for scalars
v_scalar_grid = zeros(Nsubgrid, Nsubgrid)
mu_grid = zeros(Nsubgrid, Nsubgrid)
earth_ang_sep = zeros(Nsubgrid, Nsubgrid)
moon_ang_sep = zeros(Nsubgrid, Nsubgrid)
projected_velocities_no_cb = zeros(Nsubgrid, Nsubgrid)
scattering_sub_lambert = zeros(Nsubgrid, Nsubgrid)
scattering_sub_super_lambert = zeros(Nsubgrid, Nsubgrid)

# allocate memory for RM curves
final_mean_intensity = zeros(Nt)
final_weight_v_no_cb = zeros(Nt)
final_weight_v_no_cb_lambert = zeros(Nt)
final_weight_v_no_cb_super_lambert = zeros(Nt)

# set plotting boolean
plot = false

# loop over time
for (t, epoch) in enumerate(obstimes)
    # designated time
    # if t != 8
    #     continue
    # end

    # re-zero quantities
    v_scalar_grid .= 0.0
    mu_grid .= 0.0
    earth_ang_sep .= 0.0
    moon_ang_sep .= 0.0
    projected_velocities_no_cb .= 0.0

    # get state vectors for bodies
    BE_bary = spkssb(399, epoch, "J2000")

    # get observatory position in barycenter frame
    EO_bary = sxform("ITRF93", "J2000", epoch) * EO_earth

    # get vector from barycenter to observatory on Earth's surface
    BO_bary = BE_bary .+ EO_bary

    # get light travel time corrected observer to Europa vector
    OJ2_bary, OJ2_lt, OJ2_dlt = spkltc(502, epoch, "J2000", lt_flag, BO_bary)

    # get state of europa
    BJ2_bary = spkssb(502, epoch - OJ2_lt, "J2000")

    # get state vector of bodies as seen from Europa
    J2E_bary, J2E_lt, J2E_dlt = spkltc(399, epoch - OJ2_lt, "J2000", lt_flag, BJ2_bary)
    J2M_bary, J2M_lt, J2M_dlt = spkltc(301, epoch - OJ2_lt, "J2000", lt_flag, BJ2_bary)
    J2S_bary, J2S_lt, J2S_dlt = spkltc(10, epoch - OJ2_lt, "J2000", lt_flag, BJ2_bary)

    # get rotation matrix for Sun
    sun_rot_mat = pxform("IAU_SUN", "J2000", epoch - OJ2_lt - J2S_lt)

    ###### BEGIN PLOTTING BLOCK ######
    if plot
        # set up figure
        # local fig = plt.figure()
        # local ax1 = fig.add_subplot(projection="3d")

        # set up figure
        local fig = plt.figure()
        local ax1 = fig.add_subplot()
    end
    ###### END PLOTTING BLOCK ######

    # loop over disk positions
    for i in eachindex(disk.ϕc)
        for j in 1:disk.Nθ[i]
            # save the tile position
            ϕc[i,j] = disk.ϕc[i]
            θc[i,j] = disk.θc[i,j]

            # subdivide the tile
            ϕe_sub = range(disk.ϕe[i], disk.ϕe[i+1], length=Nsubgrid+1)
            θe_sub = range(disk.θe[i,j], disk.θe[i,j+1], length=Nsubgrid+1)
            ϕc_sub = GRASS.get_grid_centers(ϕe_sub)
            θc_sub = GRASS.get_grid_centers(θe_sub)
            subgrid = Iterators.product(ϕc_sub, θc_sub)

            # get sun to patch vectors
            SP_sun_pos .= map(x -> pgrrec("SUN", getindex(x,2), getindex(x,1), 0.0, sun_radius, 0.0), subgrid)

            # get differential rotation velocities
            v_scalar_grid .= map(x -> GRASS.v_scalar(x...), subgrid)

            # convert v_scalar to from km/day km/s
            v_scalar_grid ./= 86400.0

            # determine pole vector for each patch
            GRASS.pole_vector_grid!(SP_sun_pos, pole_vector_grid)

            # get velocity vector direction and set magnitude
            GRASS.v_vector(SP_sun_pos, pole_vector_grid, v_scalar_grid, SP_sun_vel)

            # rotate sun coordinates into barycenter frame
            for k in eachindex(SP_sun_pos)
                SP_bary_pos[k] .= (sun_rot_mat * SP_sun_pos[k])
                SP_bary_vel[k] .= (sun_rot_mat * SP_sun_vel[k])
                SP_bary[k] = vcat(SP_bary_pos[k], SP_bary_vel[k])
            end

            # get vector from Europa to each patch on Sun's surface
            for k in eachindex(J2P_bary)
                J2P_bary[k] = J2S_bary .+ SP_bary[k]
            end

            # calculate mu at each point
            GRASS.calc_mu_grid!(SP_bary, J2P_bary, mu_grid)

            # move on if everything is off the grid
            all(mu_grid .< zero(eltype(mu_grid))) && continue

            # get projected velocity for each patch
            GRASS.projected!(SP_bary, J2P_bary, projected_velocities_no_cb)

            # convert from km/s to m/s
            projected_velocities_no_cb .*= 1000.0
            z_rot_sub = projected_velocities_no_cb ./ GRASS.c_ms

            # determine patches that are blocked by Earth
            # calculate the distance between tile corner and moon
            for k in eachindex(J2P_bary)
                earth_ang_sep[k] = GRASS.calc_proj_dist(J2E_bary[1:3], J2P_bary[k][1:3])
                moon_ang_sep[k] = GRASS.calc_proj_dist(J2M_bary[1:3], J2P_bary[k][1:3])
            end

            # get the scattering
            scattering_sub_lambert = lambertian.(cos.(earth_ang_sep))
            scattering_sub_super_lambert = lambertian.(cos.(earth_ang_sep)).^1e8

            # get indices for visible patches
            idx1 = mu_grid .> 0.0
            # idx3 = (idx1) .& (earth_ang_sep .> atan((earth_radius)/norm(J2E_bary[1:3]))) .& (moon_ang_sep .> atan((moon_radius)/norm(J2M_bary[1:3])))
            idx3 = (idx1) .& (earth_ang_sep .> atan((earth_radius)/norm(J2E_bary[1:3])))

            # assign the mean mu as the mean of visible mus
            μs[i,j] = mean(view(mu_grid, idx1))

            # calculate projected area element of tile
            dϕ = step(ϕe_sub)
            dθ = step(θe_sub)
            dA_sub = map(x -> GRASS.calc_dA(sun_radius, getindex(x,1), dϕ, dθ), subgrid)
            dA_sub .*= mu_grid
            dA[i,j] = sum(view(dA_sub, idx1))

            # calc limb darkening
            ld_sub = map(x -> GRASS.quad_limb_darkening_eclipse(x), mu_grid)
            ld[i,j] = mean(view(ld_sub, idx3))

            # assign the scattering as the weighted mean
            scattering_lambert[i,j] = sum(view(scattering_sub_lambert .* ld_sub .* dA_sub, idx3)) ./ sum(view(ld_sub .* dA_sub, idx3))
            scattering_super_lambert[i,j] = sum(view(scattering_sub_super_lambert .* ld_sub .* dA_sub, idx3)) ./ sum(view(ld_sub .* dA_sub, idx3))

            # calc weights
            wts[i,j] = mean(view(ld_sub .* dA_sub, idx3))

            # calc z_rot
            z_rot[i,j] = sum(view(z_rot_sub .* ld_sub .* dA_sub, idx3)) ./ sum(view(ld_sub .* dA_sub, idx3))

            # set weighted mean velocitty
            mean_weight_v_no_cb[i,j] = mean(view(projected_velocities_no_cb, idx3))

            if isnan(ld[i,j])
                ld[i,j] = 0.0
                wts[i,j] = 0.0
                z_rot[i,j] = 0.0
                mean_weight_v_no_cb[i,j] = 0.0
                scattering_lambert[i,j] = 0.0
                scattering_super_lambert[i,j] = 0.0
            end

            ###### BEGIN PLOTTING BLOCK ######
            if plot
                # convert to ra dec
                OP_ra_dec = SPICE.recrad.([x[1:3] for x in J2P_bary])
                ra = getindex.(OP_ra_dec,2)
                dec = getindex.(OP_ra_dec,3)

                # plot stuff
                x = getindex.(J2P_bary,1)
                y = getindex.(J2P_bary,2)
                z = getindex.(J2P_bary,3)

                # set values to NaN for plotting
                # mu_grid[.!idx3] .= NaN
                # projected_velocities_no_cb[.!idx3] .= NaN

                cmap = mpl.cm.viridis.copy()
                cmap.set_bad(color="black")

                # ax1.pcolormesh(rad2deg.(ra), rad2deg.(dec), mu_grid, vmin=0.0, vmax=1.0, cmap=cmap)
                # ax1.pcolormesh(rad2deg.(ra), rad2deg.(dec), ld_sub .* idx3, vmin=0.0, vmax=1.0, cmap="afmhot")
                # ax1.pcolormesh(rad2deg.(ra), rad2deg.(dec), ld_sub .* scattering_sub .* idx3, vmin=0.0, vmax=1.0, cmap="afmhot")
                ax1.pcolormesh(rad2deg.(ra), rad2deg.(dec), scattering_sub_super_lambert .* idx3, vmin=0.0, vmax=1.0, cmap="viridis")
                # ax1.pcolormesh(rad2deg.(ra), rad2deg.(dec), rad2deg.(earth_ang_sep), vmin=0.0, vmax=0.1)
                # ax1.pcolormesh(rad2deg.(ra), rad2deg.(dec), projected_velocities_no_cb .* idx3, vmin=-2000, vmax=2000, cmap="seismic")
                # ax1.plot_surface(x, y, z, facecolors=mpl.cm.viridis(mu_grid))
            end
            ###### END PLOTTING BLOCK ######

        end
    end

    ###### BEGIN PLOTTING BLOCK ######
    if plot
        ax1.invert_xaxis()
        ax1.set_aspect("equal")
        plt.show()
    end
    ###### END PLOTTING BLOCK ######

    # do disk integration
    idx_grid = ld .> 0.0

    # get flux weight
    contrast = (ld / NaNMath.maximum(ld)).^0.1
    brightness = ld .* dA
    cheapflux = sum(view(brightness, idx_grid))

    # determine final mean intensity for disk grid
    final_mean_intensity[t] = cheapflux

    #determine final mean weighted velocity for disk grid
    final_weight_v_no_cb[t] = sum(view(contrast .* mean_weight_v_no_cb .* brightness, idx_grid)) / cheapflux
    # final_weight_v_no_cb += mean(view(mean_weight_v_earth_orb, idx_grid))

    # get flux weights with scattering
    contrast = (ld / NaNMath.maximum(ld)).^0.1
    brightness = ld .* dA .* scattering
    cheapflux = sum(view(brightness, idx_grid))

    final_weight_v_no_cb_lambert[t] = sum(view(contrast .* mean_weight_v_no_cb .* brightness, idx_grid)) / cheapflux
    # final_weight_v_no_cb_lambert += mean(view(mean_weight_v_earth_orb, idx_grid))

    # get flux weights with scattering
    contrast = (ld / NaNMath.maximum(ld)).^0.1
    brightness = ld .* dA .* scattering
    cheapflux = sum(view(brightness, idx_grid))

    final_weight_v_no_cb_super_lambert[t] = sum(view(contrast .* mean_weight_v_no_cb .* brightness, idx_grid)) / cheapflux
    # final_weight_v_no_cb_super_lambert += mean(view(mean_weight_v_earth_orb, idx_grid))
end

println(">>> Plotting RM...")

# get time axis
taxis = map(x -> et2utc(x, "ISOD", 0), obstimes)

# initialize plot object
fig, ax1 = plt.subplots()

# plot time series
ax1.scatter(eachindex(taxis), final_weight_v_no_cb, label=L"{\rm Classical\ RM}")
ax2.scatter(eachindex(taxis), final_weight_v_no_cb_lambert, label=L"{\rm RM\ + Lambertian\ Reflector}")
ax2.scatter(eachindex(taxis), final_weight_v_no_cb_super_lambert, label=L"{\rm RM\ + Super\ Lambertian\ Reflector}")

# set labels and other pretty stuff
ax1.set_xlabel(L"{\rm Time}")
ax1.set_ylabel(L"{\rm RV\ (m/s)}")
ax1.legend()

plt.show()
