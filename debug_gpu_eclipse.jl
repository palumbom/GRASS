using Revise
using GRASS
using SPICE
using CSV
using DataFrames
import PyPlot
using Statistics
plt = PyPlot

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938

#convert from utc to et as needed by SPICE
time_stamps = utc2et.(["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"])

N = 50
Nt = length(time_stamps)
Nsubgrid = 8
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
gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, 1)

wsp = GRASS.SynthWorkspaceEclipse(disk, 1, Nt, verbose=true)
mem = GRASS.GeoWorkspaceEclipse(disk, 1, Nt)
LD_type = "SSD"

mu_grid_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))    
dA_total_proj_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))
idx1_matrix = Matrix{Matrix{Int}}(undef, length(disk.ϕc), maximum(disk.Nθ))
idx3_matrix = Matrix{Matrix{Int}}(undef, length(disk.ϕc), maximum(disk.Nθ))
z_rot_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))
zenith_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))

# compute geometry for timestamp
for t in 1:length(time_stamps)
    
    GRASS.eclipse_compute_quantities!(disk, time_stamps[t], t, obs_long, obs_lat, alt, wsp.ϕc, wsp.θc, 
                                    wsp.μs, wsp.z_rot, wsp.dA, wsp.xyz, wsp.ax_codes, 
                                    mem.dA_total_proj_mean, mem.mean_weight_v_no_cb,
                                    mem.mean_weight_v_earth_orb, mem.pole_vector_grid,
                                    mem.SP_sun_pos, mem.SP_sun_vel, mem.SP_bary, mem.SP_bary_pos,
                                    mem.SP_bary_vel, mem.OP_bary, mem.mu_grid, mem.projected_velocities_no_cb, 
                                    mem.distance, mem.v_scalar_grid, mem.v_earth_orb_proj, mu_grid_matrix,
                                    dA_total_proj_matrix, idx1_matrix, idx3_matrix, z_rot_matrix, zenith_matrix)    

# lambdas_cpu, outspec_cpu = GRASS.eclipse_compute_intensity(disk, lines, 1.0, LD_type, idx1_matrix, idx3_matrix,
#                                 mu_grid_matrix, z_rot_matrix, dA_total_proj_matrix, wsp.ld, wsp.z_rot, zenith_matrix, 
#                                 wsp.μs, wsp.ax_codes, wsp.dA, wsp.μs, wsp.ax_codes, wsp.dA, true, t, wsp.ext)

    mean_intensity = GRASS.eclipse_compute_intensity(disk, lines, 1.0, LD_type, idx1_matrix, idx3_matrix,
                                    mu_grid_matrix, mem.mean_weight_v_no_cb[:, :, t], mem.mean_weight_v_earth_orb[:, :, t],
                                    z_rot_matrix, dA_total_proj_matrix, wsp.ld, wsp.z_rot, zenith_matrix, true, wsp.ext)
    idx_grid = mean_intensity[:, :, 1] .> 0.0
    brightness = mean_intensity[:, :, 1] .* mem.dA_total_proj_mean[:, :, t] .* wsp.ext[:, :, 1]
    cheapflux = sum(view(brightness, idx_grid))

    final_weight_v_no_cb = sum(view(mem.mean_weight_v_no_cb[:, :, t] .* brightness, idx_grid)) / cheapflux 
    final_weight_v_no_cb += mean(view(mem.mean_weight_v_earth_orb[:, :, t], idx_grid))                       
    
    """
    GPU GEOMETRY CODE CALLS BELOW
    """

    GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, lines, LD_type, 1.0, 1.0, disk, gpu_allocs)

    idx_grid2 = Array(gpu_allocs.ld[:, :, 1]) .> 0.0

    brightness2 = Array(gpu_allocs.ld[:, :, 1]) .* Array(gpu_allocs.dA) .* Array(gpu_allocs.ext[:, :, 1])

    cheapflux2 = sum(view(brightness2, idx_grid2))

    final_weight_v_no_cb2 = sum(view(Array(gpu_allocs.projected_v) .* brightness2, idx_grid2)) / cheapflux2
    final_weight_v_no_cb2 += mean(view(Array(gpu_allocs.earth_v), idx_grid2))  

    # plt.imshow(mem.mean_weight_v_earth_orb[:, :, t] .- Array(gpu_allocs.earth_v))
    println(final_weight_v_no_cb - final_weight_v_no_cb2)
end
# plt.savefig("test")
# plt.clf()

# make a plot
# plt.imshow(wsp.z_rot .- Array(gpu_allocs.z_rot))
# plt.imshow(wsp.z_rot .* GRASS.c_ms, vmin=-2000, vmax=2000, cmap="seismic")
# plt.imshow(Array(gpu_allocs.z_rot) .* GRASS.c_ms .- wsp.z_rot .* GRASS.c_ms)
# plt.title("z_rot")
# plt.colorbar()
# plt.savefig("gpu_test.png")
# plt.show()
# plt.clf()
