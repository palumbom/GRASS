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
time_stamps = utc2et.(["2023-10-14T14:29:23.500000"])

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
t = 1
    
GRASS.eclipse_compute_quantities!(disk, time_stamps[t], t, obs_long, obs_lat, alt, wsp.ϕc, wsp.θc, 
                                  wsp.μs, wsp.z_rot, wsp.dA, wsp.xyz, wsp.ax_codes, 
                                  mem.dA_total_proj_mean, mem.mean_weight_v_no_cb,
                                  mem.mean_weight_v_earth_orb, mem.pole_vector_grid,
                                  mem.SP_sun_pos, mem.SP_sun_vel, mem.SP_bary, mem.SP_bary_pos,
                                  mem.SP_bary_vel, mem.OP_bary, mem.mu_grid, mem.projected_velocities_no_cb, 
                                  mem.distance, mem.v_scalar_grid, mem.v_earth_orb_proj, mu_grid_matrix,
                                  dA_total_proj_matrix, idx1_matrix, idx3_matrix, z_rot_matrix, zenith_matrix)    

lambdas_cpu, outspec_cpu = GRASS.eclipse_compute_intensity(disk, lines, 1.0, LD_type, idx1_matrix, idx3_matrix,
                                mu_grid_matrix, z_rot_matrix, dA_total_proj_matrix, wsp.ld, wsp.z_rot, zenith_matrix, 
                                wsp.μs, wsp.ax_codes, wsp.dA, wsp.μs, wsp.ax_codes, wsp.dA, true, t, wsp.ext)
   
"""
GPU GEOMETRY CODE CALLS BELOW
"""
gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, 1)

GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, lines, LD_type, 1.0, 1.0, disk, gpu_allocs)

# make a plot
plt.imshow(wsp.ext .- Array(gpu_allocs.ext))
# plt.imshow(wsp.z_rot .* GRASS.c_ms, vmin=-2000, vmax=2000, cmap="seismic")
# plt.imshow(Array(gpu_allocs.z_rot) .* GRASS.c_ms .- wsp.z_rot .* GRASS.c_ms)
plt.title("ext")
plt.colorbar()
plt.savefig("gpu_test.png")
plt.show()
plt.clf()
