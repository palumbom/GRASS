using Revise
using GRASS
using SPICE

# epoch = "2023-10-14T15:26:45.500000"
epoch = utc2et("2023-10-14T17:29:23.500000")

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938

disk = GRASS.DiskParamsEclipse(N=100, Nt=1, Nsubgrid=5)

lines = [15652.79] # array of line centers 
depths = [0.25]   # array of line depths
templates = ["FeI_5434"] # template data to use
variability = trues(length(lines))  # whether or not the bisectors should "dance"
blueshifts = zeros(length(lines))   # set convective blueshift value
resolution = 7e5                    # spectral resolution

# make the disk and spec composite type instances
wsp = GRASS.SynthWorkspaceEclipse(disk, 1, 1, verbose=true)
spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
                        blueshifts=blueshifts, templates=templates, resolution=resolution)  
gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, 1)

GRASS.calc_eclipse_quantities_gpu!(epoch, obs_long, obs_lat, alt, [5500.0], disk, gpu_allocs)

# get size of sub-tiled grid
Nϕ = disk.N
Nθ_max = maximum(disk.Nθ)
Nsubgrid = disk.Nsubgrid
Nϕ_sub = Nϕ * Nsubgrid
Nθ_sub = maximum(disk.Nθ) * Nsubgrid

# disk coordinates
ϕe = disk.ϕe
ϕc = disk.ϕc
θe = disk.θe
θc = disk.θc

#query JPL horizons for E, S, M position (km) and velocities (km/s)
BE_bary = spkssb(399,epoch,"J2000")

#determine xyz earth coordinates for lat/long of observatory
flat_coeff = (GRASS.earth_radius - GRASS.earth_radius_pole) / GRASS.earth_radius
EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, GRASS.earth_radius, flat_coeff)
#set earth velocity vectors
EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
#transform into ICRF frame
EO_bary = sxform("ITRF93", "J2000", epoch) * EO_earth

# get vector from barycenter to observatory on Earth's surface
BO_bary = BE_bary .+ EO_bary

# set string for ltt and abberation
lt_flag = "CN+S"

# get light travel time corrected OS vector
OS_bary, OS_lt, OS_dlt = spkltc(10, epoch, "J2000", lt_flag, BO_bary)

# get vector from observatory on earth's surface to moon center
OM_bary, OM_lt, OM_dlt = spkltc(301, epoch, "J2000", lt_flag, BO_bary)

# get modified epch
epoch_lt = epoch - OS_lt

# get rotation matrix for sun
sun_rot_mat = pxform("IAU_SUN", "J2000", epoch_lt)
R_x = sun_rot_mat


dat = Array(gpu_allocs.z_rot) .* GRASS.c_ms
# cmap = plt.cm.inferno
cmap = plt.cm.seismic
vmin = -2000.0
vmax = 2000.0

norm = plt.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

# initialize figure
fig, ax1 = plt.subplots()

# plot circle background to smooth jagged edges
circle1 = plt.matplotlib.patches.Circle((0, 0), 1.01, color="k", zorder=0)
ax1.add_patch(circle1)

# coords for zoom in cell
ϕidx = 12
θidx = 51

# loop over grid positions
println("\t>>> Plotting!")
for i in 1:length(ϕe)-1
    lat = range(ϕe[i], ϕe[i+1], length=4)
    for j in 1:disk.Nθ[i]
        lon = range(θe[i,j], θe[i,j+1], length=4)

        border = (([lat[1], lat[2], lat[3], lat[4], lat[4], lat[4], lat[4], lat[3], lat[2], lat[1], lat[1], lat[1], lat[1]]),
                  ([lon[1], lon[1], lon[1], lon[1], lon[2], lon[3], lon[4], lon[4], lon[4], lon[4], lon[3], lon[2], lon[1]]))


        out = GRASS.sphere_to_cart.(1.0, border...)
        x = getindex.(out, 1)
        y = getindex.(out, 2)
        z = getindex.(out, 3)

        # rotate it
        for k in eachindex(x)
            x0 = x[k]
            y0 = y[k]
            z0 = z[k]

            x[k] = x0 * R_x[1,1] + y0 * R_x[1,2] + z0 * R_x[1,3]
            y[k] = x0 * R_x[2,1] + y0 * R_x[2,2] + z0 * R_x[2,3]
            z[k] = x0 * R_x[3,1] + y0 * R_x[3,2] + z0 * R_x[3,3]
        end

        idx = z .>= 0

        color="k"
        lw = 1
        zorder=1

        if any(idx)
            # ax1.plot(x[idx], y[idx], color=color, lw=lw, zorder=zorder)
            ax1.fill(x[idx], y[idx], c=smap.to_rgba(dat[i,j]), zorder=0)
        end
    end
end
ax1.set_aspect("equal")
plt.show()