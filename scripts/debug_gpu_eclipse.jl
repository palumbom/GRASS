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
Nsubgrid = 10

lines = [5434.5232] # array of line centers 
depths = [0.6]   # array of line depths
templates = ["FeI_5434"] # template data to use
variability = trues(length(lines))  # whether or not the bisectors should "dance"
blueshifts = zeros(length(lines))   # set convective blueshift value
resolution = 7e5                    # spectral resolution

LD_type = "SSD_4parameter"
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=Nsubgrid)

ext_coeff_array = [0.15452995224327976, 0.15256098077094832, 0.14720055859068512, 0.154895798933504, 0.15181381895180662, 0.15107508233588227, 0.15116772762156633, 0.14882114581650618, 0.14865707189399568, 0.1494903120065096, 0.16011027092744037, 0.15593033972594958, 0.14195968590211427, 0.15401904166429853, 0.1277699772941639, 0.12709315507233226, 0.12820346527304866, 0.11702310600015708, 0.1435320747844216, 0.12380490304619193, 0.12450734135297492, 0.12101777355247835]
for t in 1:length(time_stamps)
    # make the disk and spec composite type instances
    spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability, blueshifts=blueshifts, templates=templates, resolution=resolution) 

    gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, 1)
    wsp = GRASS.SynthWorkspaceEclipse(disk, 1, Nt, verbose=true)

    """
    CPU GEOMETRY CODE CALLS BELOW
    """
    GRASS.eclipse_compute_quantities!(time_stamps[t], t, obs_long, obs_lat, alt, lines, LD_type, true, ext_coeff_array[t], disk, wsp) 

    idx_grid = wsp.ld[:, :, 1] .> 0.0
    brightness = wsp.ld[:, :, 1] .* wsp.dA[:, :, t] .* wsp.ext[:, :, 1]
    cheapflux = sum(view(brightness, idx_grid))

    final_weight_v_no_cb = sum(view(wsp.mean_weight_v_no_cb[:, :, t] .* brightness, idx_grid)) / cheapflux 
    final_weight_v_no_cb += mean(view(wsp.mean_weight_v_earth_orb[:, :, t], idx_grid))                       
    
    """
    GPU GEOMETRY CODE CALLS BELOW
    """
    GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, lines, LD_type, 1.0, ext_coeff_array[t], disk, gpu_allocs, 0.0)

    idx_grid2 = Array(gpu_allocs.ld[:, :, 1]) .> 0.0
    brightness2 = Array(gpu_allocs.ld[:, :, 1]) .* Array(gpu_allocs.dA) .* Array(gpu_allocs.ext[:, :, 1])
    cheapflux2 = sum(view(brightness2, idx_grid2))

    final_weight_v_no_cb2 = sum(view(Array(gpu_allocs.projected_v) .* brightness2, idx_grid2)) / cheapflux2
    final_weight_v_no_cb2 += mean(view(Array(gpu_allocs.earth_v), idx_grid2))  

    plt.imshow(wsp.mean_weight_v_earth_orb[:, :, t] .- Array(gpu_allocs.earth_v))
end
plt.colorbar()
plt.savefig("debug_gpu_test.png")
plt.show()
plt.clf()

# EO_bary = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# moon_radius = 1.0
# sun_radius = 1.0