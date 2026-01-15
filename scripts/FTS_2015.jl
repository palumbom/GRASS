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

GRASS.get_kernels()

# determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)

function all_lines_gpu(time, granulation_status, LD_type, ext_toggle, model, spot_toggle)
    # convert from utc to et as needed by SPICE
    time_stamps = utc2et.(time)

    # FTS location
    obs_lat = 51.560583 
    obs_long = 9.944333  
    alt = 0.201 

    # set up paramaters for disk
    N = 197
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=40)

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    depth = GRASS.get_depth(lp)
    lfile = GRASS.get_file(lp)
    #original
    rv = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error = Vector{Vector{Float64}}(undef,size(name)...)
    resolution = 7e5
    
    for i in eachindex(lp.λrest) 

        # set up parameters for synthesis
        lines = [λrest[i]]
        templates = [lfile[i]]
        depths = [depth[i]]

        if granulation_status == true
            variability = trues(length(lines))  # whether or not the bisectors should "dance"
        end
        if granulation_status == false
            variability = falses(length(lines))
        end
        blueshifts = zeros(length(lines))   # set convective blueshift value

        # make the spec composite type instances
        spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
                        blueshifts=blueshifts, templates=templates, resolution=resolution) 

        if model == "LD" 
            lambdas, outspec = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, obs_long, obs_lat, alt, time_stamps, lines, 0.0, ext_toggle, spot_toggle)
        end

        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas, outspec, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)
    end
    return rv, rv_error, λrest
end

function FTS_eclipse_var_on_gpu(LD_type, ext_toggle, model, spot_toggle)
    FTS = ["2015-03-20T7:07:57.2", "2015-03-20T7:09:45.8", "2015-03-20T7:11:34.3", "2015-03-20T7:13:22.9", "2015-03-20T7:15:11.5", "2015-03-20T7:17:00.3", "2015-03-20T7:18:49.7", "2015-03-20T7:20:38.5", "2015-03-20T7:22:27.4", "2015-03-20T7:24:16.5", "2015-03-20T7:26:05.5", "2015-03-20T7:27:53.8", "2015-03-20T7:29:42.5", "2015-03-20T7:31:30.7", "2015-03-20T7:33:19.2", "2015-03-20T7:35:09.5", "2015-03-20T7:36:58.0", "2015-03-20T7:38:46.1", "2015-03-20T7:40:34.6", "2015-03-20T7:42:22.8", "2015-03-20T7:44:11.5", "2015-03-20T7:46:00.1", "2015-03-20T7:47:48.8", "2015-03-20T7:49:37.1", "2015-03-20T7:51:25.3", "2015-03-20T7:53:13.9", "2015-03-20T7:55:02.3", "2015-03-20T7:56:50.7", "2015-03-20T7:58:39.0", "2015-03-20T8:00:27.2", "2015-03-20T8:02:15.7", "2015-03-20T8:04:04.0", "2015-03-20T8:05:53.0", "2015-03-20T8:07:41.4", "2015-03-20T8:09:30.0", "2015-03-20T8:11:18.5", "2015-03-20T8:13:07.1", "2015-03-20T8:14:55.3", "2015-03-20T8:16:43.7", "2015-03-20T8:18:32.1", "2015-03-20T8:20:20.4", "2015-03-20T8:22:09.4", "2015-03-20T8:23:57.7", "2015-03-20T8:25:46.1", "2015-03-20T8:27:34.7", "2015-03-20T8:29:23.0", "2015-03-20T8:31:11.3", "2015-03-20T8:32:59.5", "2015-03-20T8:34:47.9", "2015-03-20T8:36:36.3", "2015-03-20T8:38:54.5", "2015-03-20T8:40:43.0", "2015-03-20T8:42:31.4", "2015-03-20T8:44:19.9", "2015-03-20T8:46:08.5", "2015-03-20T8:47:56.8", "2015-03-20T8:49:45.5", "2015-03-20T8:51:34.0", "2015-03-20T8:53:22.9", "2015-03-20T8:55:11.5", "2015-03-20T8:56:59.6", "2015-03-20T8:58:47.9", "2015-03-20T9:00:36.1", "2015-03-20T9:02:24.8", "2015-03-20T9:04:13.1", "2015-03-20T9:06:01.7", "2015-03-20T9:07:50.2", "2015-03-20T9:09:38.7", "2015-03-20T9:11:27.1", "2015-03-20T9:13:15.5", "2015-03-20T9:15:04.0", "2015-03-20T9:16:53.1", "2015-03-20T9:18:42.1", "2015-03-20T9:20:30.9", "2015-03-20T9:22:19.7", "2015-03-20T9:24:08.5", "2015-03-20T9:25:57.3", "2015-03-20T9:27:45.9", "2015-03-20T9:29:34.7", "2015-03-20T9:31:23.1", "2015-03-20T9:33:11.5", "2015-03-20T9:34:59.9", "2015-03-20T9:36:48.4", "2015-03-20T9:38:37.4", "2015-03-20T9:40:26.3", "2015-03-20T9:42:15.0", "2015-03-20T9:44:03.8", "2015-03-20T9:45:52.2", "2015-03-20T9:47:40.7", "2015-03-20T9:49:29.7", "2015-03-20T9:51:18.4", "2015-03-20T9:53:07.5", "2015-03-20T9:54:55.9", "2015-03-20T9:56:44.9", "2015-03-20T9:58:33.3", "2015-03-20T10:00:21.8", "2015-03-20T10:02:10.1", "2015-03-20T10:03:58.6", "2015-03-20T10:05:47.4", "2015-03-20T10:07:36.2", "2015-03-20T10:09:54.5", "2015-03-20T10:11:43.9", "2015-03-20T10:13:33.6", "2015-03-20T10:15:22.6", "2015-03-20T10:17:11.7", "2015-03-20T10:19:00.9", "2015-03-20T10:20:49.9", "2015-03-20T10:22:38.9", "2015-03-20T10:24:27.9", "2015-03-20T10:26:17.0", "2015-03-20T10:28:07.1", "2015-03-20T10:29:56.1", "2015-03-20T10:31:45.1", "2015-03-20T10:33:34.0", "2015-03-20T10:35:22.9", "2015-03-20T10:37:12.0", "2015-03-20T10:39:01.0", "2015-03-20T10:40:49.9", "2015-03-20T10:42:38.7", "2015-03-20T10:44:27.6", "2015-03-20T10:46:16.8", "2015-03-20T10:48:05.8", "2015-03-20T10:49:54.9", "2015-03-20T10:51:43.6", "2015-03-20T10:53:32.7", "2015-03-20T10:55:21.9", "2015-03-20T10:57:10.8", "2015-03-20T10:58:59.7", "2015-03-20T11:00:49.8", "2015-03-20T11:02:38.6", "2015-03-20T11:04:27.7", "2015-03-20T11:06:16.8", "2015-03-20T11:08:05.7", "2015-03-20T11:09:54.6", "2015-03-20T11:11:43.6", "2015-03-20T11:13:33.1", "2015-03-20T11:15:22.0", "2015-03-20T11:17:10.9", "2015-03-20T11:18:59.9", "2015-03-20T11:20:48.7", "2015-03-20T11:22:37.7", "2015-03-20T11:24:26.7", "2015-03-20T11:26:15.7", "2015-03-20T11:28:04.3", "2015-03-20T11:29:53.0", "2015-03-20T11:31:41.7", "2015-03-20T11:33:30.3", "2015-03-20T11:35:19.0", "2015-03-20T11:37:07.7", "2015-03-20T11:38:56.9", "2015-03-20T11:48:37.5", "2015-03-20T11:50:26.6", "2015-03-20T11:52:15.6", "2015-03-20T11:54:04.2", "2015-03-20T11:55:53.1", "2015-03-20T11:57:41.8", "2015-03-20T11:59:30.5", "2015-03-20T12:01:19.1", "2015-03-20T12:03:07.6"]

    rv, rv_error, λrest = deepcopy(all_lines_gpu(FTS, true, LD_type, ext_toggle, model, spot_toggle))
    if model == "LD" 
        @save "all_lines_rv_on_$(LD_type)_gpu.jld2"
        jldopen("all_lines_rv_on_$(LD_type)_gpu.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
        end
    end
end

FTS_eclipse_var_on_gpu("SSD_4parameter", false, "LD", false)