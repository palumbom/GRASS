using CSV
using CUDA
using FFTW
using JLD2
using GRASS
using FileIO
using Printf
using Revise
using LsqFit
using SPICE
using DataFrames
using Statistics
using EchelleCCFs
using Polynomials
using BenchmarkTools
using HypothesisTests
using Interpolations
using EchelleCCFs: λ_air_to_vac, λ_vac_to_air
using FITSIO
using RvSpectMLBase
using EchelleInstruments
using Dierckx
using Dates 
using Optimization, OptimizationOptimJL

using PyCall
py"""
import sys
sys.path.append(".")
"""
datetime = pyimport("datetime")

GRASS.get_kernels()

# convert from utc to et as needed by SPICE
neid_time = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
neid_time = neid_time[1:130]
time_stamps = utc2et.(neid_time)
neid_time = [datetime.datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f") for t in neid_time]
timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
timestamps = timestamps_full_october[16:length(timestamps_full_october)-150]
timestamps = timestamps[1:130]
path = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"

# NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

line_data = jldopen("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "r") do file
    Dict(var => read(file, var) for var in ["rv"])
end
observed_data = line_data["rv"][6] 
observed_data = observed_data .- observed_data[length(observed_data)]

rv = Vector{Float64}(undef,size(time_stamps)...)
# get lines to construct templates
lp = GRASS.LineProperties()
λrest = GRASS.get_rest_wavelength(lp)
depth = GRASS.get_depth(lp)
lfile = GRASS.get_file(lp)
resolution = 7e5

# set up paramaters for disk
N = 197
Nt = length(time_stamps)
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=40)

ext_coeff_array = [0.15452995224327976, 0.15256098077094832, 0.14720055859068512, 0.154895798933504, 0.15181381895180662, 0.15107508233588227, 0.15116772762156633, 0.14882114581650618, 0.14865707189399568, 0.1494903120065096, 0.16011027092744037, 0.15593033972594958, 0.14195968590211427, 0.15401904166429853, 0.1277699772941639, 0.12709315507233226, 0.12820346527304866, 0.11702310600015708, 0.1435320747844216, 0.12380490304619193, 0.12450734135297492, 0.12101777355247835]

function projected_RV_gpu(params)
    for i in 6:6
        println("\t>>> Template: " * string(splitdir(lfile[i])[2]))
        # set up parameters for synthesis
        lines = [λrest[i]]
        templates = [lfile[i]]
        depths = [depth[i]]
        variability = trues(length(lines)) 
        blueshifts = zeros(length(lines))   # set convective blueshift value

        # make the spec composite type instances
        spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
                        blueshifts=blueshifts, templates=templates, resolution=resolution) 
        for t in 1:disk.Nt

            gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, Int(length(lines)), precision=Float64, verbose=true)
            GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, lines, ext_coeff_array[i], disk, gpu_allocs,
                                                params[1], params[2], params[3], params[4], params[5])

            idx_grid = Array(gpu_allocs.ld[:, :, 1]) .> 0.0

            brightness = Array(gpu_allocs.ld[:, :, 1]) .* Array(gpu_allocs.dA) .* Array(gpu_allocs.ext[:, :, 1])

            cheapflux = sum(view(brightness, idx_grid))

            # determine final mean weighted velocity for disk grid
            final_weight_v_no_cb = sum(view(Array(gpu_allocs.projected_v) .* brightness, idx_grid)) / cheapflux 
            final_weight_v_no_cb += mean(view(Array(gpu_allocs.earth_v), idx_grid))  

            rv[t] = deepcopy(final_weight_v_no_cb)
        end
    end
    return rv .- rv[length(rv)]
end

function rms(params, _)
    model_output = projected_RV_gpu(params)
    loss = sqrt(sum((observed_data .- model_output).^2) / length(observed_data))

    return loss 
end

x0 = [0.31030251182557, -1.4138084241294868, -28.857187676819656, 393.9453708988727, -542.4932089212795]
prob = Optimization.OptimizationProblem(rms, x0, nothing)
sol = solve(prob, Optim.NelderMead())
print(sol)