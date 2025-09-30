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

GRASS.get_kernels()

#convert from utc to et as needed by SPICE
time_stamps = utc2et.(["2023-10-14T15:26:45.500000"]) 

N = 50
Nt = length(time_stamps)
Nsubgrid = 8

lines = [5434.5232] # array of line centers 
depths = [0.6]   # array of line depths
templates = ["FeI_5434"] # template data to use
variability = trues(length(lines))  # whether or not the bisectors should "dance"
blueshifts = zeros(length(lines))   # set convective blueshift value
resolution = 7e5                    # spectral resolution

LD_type = "SSD_quadratic"
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=Nsubgrid)

mu_grid_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))    
dA_total_proj_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))
idx1_matrix = Matrix{Matrix{Int}}(undef, length(disk.ϕc), maximum(disk.Nθ))
idx3_matrix = Matrix{Matrix{Int}}(undef, length(disk.ϕc), maximum(disk.Nθ))
z_rot_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))
zenith_matrix = Matrix{Matrix{Float64}}(undef, length(disk.ϕc), maximum(disk.Nθ))
extinction_coeff = GRASS.ext_coeff_file
# compute geometry for timestamp
for t in 1:length(time_stamps)

    if t < 46
        neid_ext_coeff = extinction_coeff[extinction_coeff[!, "Wavelength"] .== 5434.5232, "ext1"][1]
    elseif t >= 46 
        neid_ext_coeff = extinction_coeff[extinction_coeff[!, "Wavelength"] .== 5434.5232, "ext2"][1]
    end

    # make the disk and spec composite type instances
    spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
    blueshifts=blueshifts, templates=templates, resolution=resolution) 

    gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, 1)

    wsp = GRASS.SynthWorkspaceEclipse(disk, 1, Nt, verbose=true)
    mem = GRASS.GeoWorkspaceEclipse(disk, 1, Nt)
    """
    CPU GEOMETRY CODE CALLS BELOW
    """
    
    GRASS.eclipse_compute_quantities!(disk, time_stamps[t], t, obs_long, obs_lat, alt, wsp.ϕc, wsp.θc, 
                                    wsp.μs, wsp.z_rot, wsp.dA, wsp.xyz, wsp.ax_codes, 
                                    mem.dA_total_proj_mean, mem.mean_weight_v_no_cb,
                                    mem.mean_weight_v_earth_orb, mem.pole_vector_grid,
                                    mem.SP_sun_pos, mem.SP_sun_vel, mem.SP_bary, mem.SP_bary_pos,
                                    mem.SP_bary_vel, mem.OP_bary, mem.mu_grid, mem.projected_velocities_no_cb, 
                                    mem.distance, mem.v_scalar_grid, mem.v_earth_orb_proj, mu_grid_matrix,
                                    dA_total_proj_matrix, idx1_matrix, idx3_matrix, z_rot_matrix, zenith_matrix)    

    mean_intensity = GRASS.eclipse_compute_intensity(disk, lines, neid_ext_coeff, LD_type, idx1_matrix, idx3_matrix,
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

    GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, lines, LD_type, 1.0, neid_ext_coeff, disk, gpu_allocs)

    idx_grid2 = Array(gpu_allocs.ld[:, :, 1]) .> 0.0

    brightness2 = Array(gpu_allocs.ld[:, :, 1]) .* Array(gpu_allocs.dA) .* Array(gpu_allocs.ext[:, :, 1])

    cheapflux2 = sum(view(brightness2, idx_grid2))

    final_weight_v_no_cb2 = sum(view(Array(gpu_allocs.projected_v) .* brightness2, idx_grid2)) / cheapflux2
    final_weight_v_no_cb2 += mean(view(Array(gpu_allocs.earth_v), idx_grid2))  

    plt.imshow(mem.mean_weight_v_no_cb[:, :, t] .- Array(gpu_allocs.projected_v))
end

# make a plot
plt.title("projected_v")
plt.colorbar()
plt.savefig("gpu_test.png")
plt.show()
plt.clf()
