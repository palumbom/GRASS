using GRASS
using Random
using Revise
using Statistics
using SPICE
using Pkg
using PyPlot
mpl = plt.matplotlib

##TO DO: figure out what rotation matrix to use

GRASS.get_kernels()

neid_timestamps = ["2023-10-14T15:26:45", "2023-10-14T15:28:07", "2023-10-14T15:29:30", "2023-10-14T15:30:53", "2023-10-14T15:32:15", "2023-10-14T15:33:38", "2023-10-14T15:35:01", "2023-10-14T15:36:23", "2023-10-14T15:37:46", "2023-10-14T15:39:09", "2023-10-14T15:40:31", "2023-10-14T15:41:54", "2023-10-14T15:43:17", "2023-10-14T15:44:39", "2023-10-14T15:46:02", "2023-10-14T15:47:25", "2023-10-14T15:48:47", "2023-10-14T15:50:10", "2023-10-14T15:51:33", "2023-10-14T15:52:56", "2023-10-14T15:54:18", "2023-10-14T15:55:41", "2023-10-14T15:57:04", "2023-10-14T15:58:26", "2023-10-14T15:59:49", "2023-10-14T16:01:12", "2023-10-14T16:02:34", "2023-10-14T16:03:57", "2023-10-14T16:05:20", "2023-10-14T16:06:42", "2023-10-14T16:08:05", "2023-10-14T16:09:28", "2023-10-14T16:10:50", "2023-10-14T16:12:13", "2023-10-14T16:13:36", "2023-10-14T16:14:58", "2023-10-14T16:16:21", "2023-10-14T16:17:44", "2023-10-14T16:19:06", "2023-10-14T16:20:29", "2023-10-14T16:21:52", "2023-10-14T16:23:15", "2023-10-14T16:24:37", "2023-10-14T16:26:00", "2023-10-14T16:27:23", "2023-10-14T16:28:45", "2023-10-14T16:30:08", "2023-10-14T16:31:31", "2023-10-14T16:32:53", "2023-10-14T16:34:16", "2023-10-14T16:35:39", "2023-10-14T16:37:01", "2023-10-14T16:38:24", "2023-10-14T16:39:47", "2023-10-14T16:41:09", "2023-10-14T16:42:32", "2023-10-14T16:43:55", "2023-10-14T16:45:17", "2023-10-14T16:46:40", "2023-10-14T16:48:03", "2023-10-14T16:49:25", "2023-10-14T16:50:48", "2023-10-14T16:52:11", "2023-10-14T16:53:33", "2023-10-14T16:54:56", "2023-10-14T16:56:19", "2023-10-14T16:57:42", "2023-10-14T16:59:04", "2023-10-14T17:00:27", "2023-10-14T17:01:50", "2023-10-14T17:03:12", "2023-10-14T17:04:35", "2023-10-14T17:05:58", "2023-10-14T17:07:20", "2023-10-14T17:08:43", "2023-10-14T17:10:06", "2023-10-14T17:11:28", "2023-10-14T17:12:51", "2023-10-14T17:14:14", "2023-10-14T17:15:36", "2023-10-14T17:16:59", "2023-10-14T17:18:22", "2023-10-14T17:19:44", "2023-10-14T17:21:07", "2023-10-14T17:22:30", "2023-10-14T17:23:52", "2023-10-14T17:25:15", "2023-10-14T17:26:38", "2023-10-14T17:28:01", "2023-10-14T17:29:23", "2023-10-14T17:30:46", "2023-10-14T17:32:09", "2023-10-14T17:33:31", "2023-10-14T17:34:54", "2023-10-14T17:36:17", "2023-10-14T17:37:39", "2023-10-14T17:39:02", "2023-10-14T17:40:25", "2023-10-14T17:41:47", "2023-10-14T17:43:10", "2023-10-14T17:44:33", "2023-10-14T17:45:55", "2023-10-14T17:47:18", "2023-10-14T17:48:41", "2023-10-14T17:50:03", "2023-10-14T17:51:26", "2023-10-14T17:52:49", "2023-10-14T17:54:11", "2023-10-14T17:55:34", "2023-10-14T17:56:57", "2023-10-14T17:58:20", "2023-10-14T17:59:42", "2023-10-14T18:01:05", "2023-10-14T18:02:28", "2023-10-14T18:03:50", "2023-10-14T18:05:13", "2023-10-14T18:06:36", "2023-10-14T18:07:58", "2023-10-14T18:09:21", "2023-10-14T18:10:44", "2023-10-14T18:12:06", "2023-10-14T18:13:29", "2023-10-14T18:14:52", "2023-10-14T18:16:14", "2023-10-14T18:17:37", "2023-10-14T18:19:00", "2023-10-14T18:20:22", "2023-10-14T18:21:45", "2023-10-14T18:23:08", "2023-10-14T18:24:30", "2023-10-14T18:25:53", "2023-10-14T18:27:16", "2023-10-14T18:28:38", "2023-10-14T18:30:01", "2023-10-14T18:31:24", "2023-10-14T18:32:47", "2023-10-14T18:34:09", "2023-10-14T18:35:32", "2023-10-14T18:36:55", "2023-10-14T18:38:17", "2023-10-14T18:39:40", "2023-10-14T18:41:03", "2023-10-14T18:42:25", "2023-10-14T18:43:48", "2023-10-14T18:45:11", "2023-10-14T18:46:33", "2023-10-14T18:47:56", "2023-10-14T18:49:19", "2023-10-14T18:50:41", "2023-10-14T18:52:04", "2023-10-14T18:53:27", "2023-10-14T18:54:49", "2023-10-14T18:56:12", "2023-10-14T18:57:35", "2023-10-14T18:58:57", "2023-10-14T19:00:20", "2023-10-14T19:01:43", "2023-10-14T19:03:06"]
#convert from utc to et as needed by SPICE
time_stamps = utc2et.(neid_timestamps)

# set up paramaters for spectrum
N = 50
Nt = length(time_stamps)
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

#NEID location 
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

sun_rot_mat = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
# allocate the memory for keys, velocities, ld, etc.
ϕc = zeros(size(disk.θc))
θc = zeros(size(disk.θc))
μs = zeros(size(disk.θc))
ld = zeros(size(disk.θc))
dA = zeros(size(disk.θc))
xyz = zeros(size(disk.θc)..., 3)
wts = zeros(size(disk.θc))
z_rot = zeros(size(disk.θc))
ax_codes = zeros(Int, size(disk.θc))
for i in 1:1#length(time_stamps)
    # precompute quantities
    sun_rot_mat[i] = GRASS.eclipse_compute_quantities!(disk, time_stamps[i], obs_long, obs_lat, alt, ϕc, θc, μs, ld, dA, xyz, wts, z_rot, ax_codes)
end

# disk coordinates
ϕe = disk.ϕe
ϕc = disk.ϕc
θe = disk.θe
θc = disk.θc
R_x = sun_rot_mat[1] 

# get color scalar mappable
dat = wts ./ maximum(wts)
cmap = plt.cm.afmhot

dat = z_rot .* 3e8
cmap = plt.cm.seismic

norm = mpl.colors.Normalize(vmin=minimum(dat), vmax=maximum(dat))
smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

# initialize figure
fig, ax = plt.subplots(1,1, figsize=(8,8))

# loop over grid positions
println("\t>>> Plotting!")
for i in 1:length(ϕe)-1
    lat = range(ϕe[i], ϕe[i+1], length=4)
    for j in 1:disk.Nθ[i]
        lon = range(θe[i,j], θe[i,j+1], length=4)

        border = (([lat[1], lat[2], lat[3], lat[4], lat[4], lat[4], lat[4], lat[3], lat[2], lat[1], lat[1], lat[1]]),
                  ([lon[1], lon[1], lon[1], lon[1], lon[2], lon[3], lon[4], lon[4], lon[4], lon[4], lon[3], lon[2]]))


        out = GRASS.sphere_to_cart_eclipse.(1.0, border...) 
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
        if any(idx)
            ax.fill(x[idx], y[idx], c=smap.to_rgba(dat[i,j]))
        end
    end
end

# get equator coords
latitude = deg2rad(0.0)
longitude = deg2rad.(range(0.0, 360.0, length=200))
x_eq = []
y_eq = []
z_eq = []
for i in eachindex(longitude)
    out = GRASS.sphere_to_cart_eclipse.(1.0, latitude, longitude[i])
    x = getindex(out, 1)
    y = getindex(out, 2)
    z = getindex(out, 3)

    x0 = x
    y0 = y
    z0 = z

    x = x0 * R_x[1,1] + y0 * R_x[1,2] + z0 * R_x[1,3]
    y = x0 * R_x[2,1] + y0 * R_x[2,2] + z0 * R_x[2,3]
    z = x0 * R_x[3,1] + y0 * R_x[3,2] + z0 * R_x[3,3]

    push!(x_eq, x)
    push!(y_eq, y)
    push!(z_eq, z)
end

# sort the values on increasing x
idx_eq = sortperm(x_eq)
x_eq = x_eq[idx_eq]
y_eq = y_eq[idx_eq]
z_eq = z_eq[idx_eq]

idx_eq = z_eq .> 0.0

# plot the equator
ax.plot(x_eq[idx_eq], y_eq[idx_eq], color="white", ls="--", zorder=3, alpha=0.75)

# get meridians
latitude = deg2rad.(range(-89.0, 89.0, length=200))
longitude = deg2rad.(range(0.0, 360.0, step=90.0))

for j in eachindex(longitude)
    out = GRASS.sphere_to_cart_eclipse.(1.0, latitude, longitude[j])

    out = hcat(out...)

    x = out[1,:]
    y = out[2,:]
    z = out[3,:]

    x0 = x
    y0 = y
    z0 = z

    x = x0 .* R_x[1,1] .+ y0 .* R_x[1,2] .+ z0 .* R_x[1,3]
    y = x0 .* R_x[2,1] .+ y0 .* R_x[2,2] .+ z0 .* R_x[2,3]
    z = x0 .* R_x[3,1] .+ y0 .* R_x[3,2] .+ z0 .* R_x[3,3]

    # plot the meridian
    idx = z .> 0.0
    ax.plot(x[idx], y[idx], color="white", ls="--", zorder=3, alpha=0.75)
end

ax.set_xlabel(L"\Delta {\rm x\ [Stellar\ Radii]}")
ax.set_ylabel(L"\Delta {\rm y\ [Stellar\ Radii]}")
ax.set_aspect("equal")
ax.grid(false)
ax.invert_xaxis()
cb = fig.colorbar(smap, ax=ax, fraction=0.1, shrink=0.8)
cb.set_label(L"{\rm Weighted\ Relative\ Intensity}")
plt.show()