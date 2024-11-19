using Revise
using GRASS
using SPICE
using CSV
using DataFrames
import PyPlot
plt = PyPlot

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938

#convert from utc to et as needed by SPICE
time_stamps = utc2et.([""2023-10-14T14:29:23.500000""])
# time_stamps = utc2et.(["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"])

N = 50
Nt = length(time_stamps)
Nsubgrid = 4
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=Nsubgrid)

lines = [5434.5232] # array of line centers 
depths = [0.6]   # array of line depths
templates = ["FeI_5434"] # template data to use
variability = trues(length(lines))  # whether or not the bisectors should "dance"
blueshifts = zeros(length(lines))   # set convective blueshift value
resolution = 7e5                    # spectral resolution

# make the disk and spec composite type instances
spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
                        blueshifts=blueshifts, templates=templates, resolution=resolution) 

"""
CPU GEOMETRY CODE CALLS BELOW
"""

wsp = GRASS.SynthWorkspaceEclipse(disk, 1, Nt, verbose=true)
mem = GRASS.GeoWorkspaceEclipse(disk, 1, Nt)
LD_type = "KSSD"

mu_grid_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))    
dA_total_proj_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))
idx1_matrix = Matrix{Matrix{Int}}(undef, length(disk.ϕc), maximum(disk.Nθ))
idx3_matrix = Matrix{Matrix{Int}}(undef, length(disk.ϕc), maximum(disk.Nθ))
z_rot_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))
zenith_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))

# # compute geometry for timestamp
# t = 1
    
# GRASS.eclipse_compute_quantities!(disk, time_stamps[t], obs_long, obs_lat, alt, wsp.ϕc, wsp.θc, 
#                                   wsp.μs, wsp.dA, wsp.xyz, wsp.ax_codes, 
#                                   mem.dA_total_proj_mean, mem.mean_weight_v_no_cb,
#                                   mem.mean_weight_v_earth_orb, mem.pole_vector_grid,
#                                   mem.SP_sun_pos, mem.SP_sun_vel, mem.SP_bary, mem.SP_bary_pos,
#                                   mem.SP_bary_vel, mem.OP_bary, mem.mu_grid, mem.projected_velocities_no_cb, 
#                                   mem.distance, mem.v_scalar_grid, mem.v_earth_orb_proj, t, mu_grid_matrix,
#                                   dA_total_proj_matrix, idx1_matrix, idx3_matrix, z_rot_matrix, zenith_matrix)           

# # extinction_coeff = DataFrame(CSV.File("data/NEID_three_extinction.csv"))
# # neid_ext_coeff = extinction_coeff[extinction_coeff[!, "Wavelength"] .== lines[1], "Ext1"]
# lambdas_cpu, outspec_cpu = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_type, obs_long, obs_lat, alt, time_stamps,
#                                 [zenith_matrix], [dA_total_proj_matrix], [idx1_matrix], [idx3_matrix], [mu_grid_matrix], [z_rot_matrix],
#                                 wsp.μs, wsp.ax_codes, wsp.dA, "three", ext_toggle = false, verbose=true, use_gpu=false)
   
"""
GPU GEOMETRY CODE CALLS BELOW
"""

lambdas_gpu, outspec_gpu = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_type, obs_long, obs_lat, alt, time_stamps,
                                [zenith_matrix], [dA_total_proj_matrix], [idx1_matrix], [idx3_matrix], [mu_grid_matrix], [z_rot_matrix],
                                wsp.μs, wsp.ax_codes, wsp.dA, "three", ext_toggle = false, verbose=true, use_gpu=true)
# v_grid_gpu, ccf_gpu = (GRASS.calc_ccf(lambdas_gpu, outspec_gpu, lines, depths, 11e4))
# rvs_gpu, sigs_gpu = (GRASS.calc_rvs_from_ccf(v_grid_gpu, ccf_gpu)) 

# # plt.plot(lambdas_cpu, outspec_cpu, label = "cpu")
# # plt.plot(lambdas_gpu, outspec_gpu, label = "gpu")
# plt.plot(time_stamps, rvs_gpu, label = "gpu")
# plt.legend()
# plt.savefig("gpu_test.png")
# plt.clf()


# gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, 1)

# GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, lines, LD_type, 1.0, neid_ext_coeff[1], disk, gpu_allocs)

# # make a plot
# plt.imshow(wsp.z_rot .- Array(gpu_allocs.z_rot))
# # plt.imshow(wsp.z_rot .* GRASS.c_ms, vmin=-2000, vmax=2000, cmap="seismic")
# # plt.imshow(Array(gpu_allocs.z_rot) .* GRASS.c_ms .- wsp.z_rot .* GRASS.c_ms)
# plt.colorbar()
# plt.savefig("gpu_test.png")
# plt.clf()
# #z_rot & dA


# # get size of sub-tiled grid
# Nϕ = disk.N
# Nθ_max = maximum(disk.Nθ)
# Nsubgrid = disk.Nsubgrid
# Nϕ_sub = Nϕ * Nsubgrid
# Nθ_sub = maximum(disk.Nθ) * Nsubgrid

# # disk coordinates
# ϕe = disk.ϕe
# ϕc = disk.ϕc
# θe = disk.θe
# θc = disk.θc

# #query JPL horizons for E, S, M position (km) and velocities (km/s)
# BE_bary = spkssb(399,epoch,"J2000")

# #determine xyz earth coordinates for lat/long of observatory
# flat_coeff = (GRASS.earth_radius - GRASS.earth_radius_pole) / GRASS.earth_radius
# EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, GRASS.earth_radius, flat_coeff)
# #set earth velocity vectors
# EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
# #transform into ICRF frame
# EO_bary = sxform("ITRF93", "J2000", epoch) * EO_earth

# # get vector from barycenter to observatory on Earth's surface
# BO_bary = BE_bary .+ EO_bary

# # set string for ltt and abberation
# lt_flag = "CN+S"

# # get light travel time corrected OS vector
# OS_bary, OS_lt, OS_dlt = spkltc(10, epoch, "J2000", lt_flag, BO_bary)

# # get vector from observatory on earth's surface to moon center
# OM_bary, OM_lt, OM_dlt = spkltc(301, epoch, "J2000", lt_flag, BO_bary)

# # get modified epch
# epoch_lt = epoch - OS_lt

# # get rotation matrix for sun
# sun_rot_mat = pxform("IAU_SUN", "J2000", epoch_lt)
# R_x = sun_rot_mat


# dat = Array(gpu_allocs.z_rot) .* GRASS.c_ms
# # cmap = plt.cm.inferno
# cmap = plt.cm.seismic
# vmin = -2000.0
# vmax = 2000.0

# norm = plt.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
# smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

# # initialize figure
# fig, ax1 = plt.subplots()

# # plot circle background to smooth jagged edges
# circle1 = plt.matplotlib.patches.Circle((0, 0), 1.01, color="k", zorder=0)
# ax1.add_patch(circle1)

# # coords for zoom in cell
# ϕidx = 12
# θidx = 51

# # loop over grid positions
# println("\t>>> Plotting!")
# for i in 1:length(ϕe)-1
#     lat = range(ϕe[i], ϕe[i+1], length=4)
#     for j in 1:disk.Nθ[i]
#         lon = range(θe[i,j], θe[i,j+1], length=4)

#         border = (([lat[1], lat[2], lat[3], lat[4], lat[4], lat[4], lat[4], lat[3], lat[2], lat[1], lat[1], lat[1], lat[1]]),
#                   ([lon[1], lon[1], lon[1], lon[1], lon[2], lon[3], lon[4], lon[4], lon[4], lon[4], lon[3], lon[2], lon[1]]))


#         out = GRASS.sphere_to_cart.(1.0, border...)
#         x = getindex.(out, 1)
#         y = getindex.(out, 2)
#         z = getindex.(out, 3)

#         # rotate it
#         for k in eachindex(x)
#             x0 = x[k]
#             y0 = y[k]
#             z0 = z[k]

#             x[k] = x0 * R_x[1,1] + y0 * R_x[1,2] + z0 * R_x[1,3]
#             y[k] = x0 * R_x[2,1] + y0 * R_x[2,2] + z0 * R_x[2,3]
#             z[k] = x0 * R_x[3,1] + y0 * R_x[3,2] + z0 * R_x[3,3]
#         end

#         idx = z .>= 0

#         color="k"
#         lw = 1
#         zorder=1

#         if any(idx)
#             # ax1.plot(x[idx], y[idx], color=color, lw=lw, zorder=zorder)
#             ax1.fill(x[idx], y[idx], c=smap.to_rgba(dat[i,j]), zorder=0)
#         end
#     end
# end
# ax1.set_aspect("equal")
# plt.show()