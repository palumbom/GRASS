using CUDA
using GRASS
using CSV
using FFTW
using JLD2
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
using EchelleCCFs: λ_air_to_vac

using PyCall
py"""
import sys
sys.path.append(".")
"""
neid_asymmetric_lsf = pyimport("NeidLsf")

variable_names = ["lambda_min", "lambda_max", "flux_min", "pixel_mean", "pixels", "neid_wavelength"]
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)

h5_path = "/storage/home/efg5335/work/sw/DISCO/GRASS/disco_fe5250_params.h5"

ext_coeff_array = [0.1438304560465685, 0.14187742832512693, 0.13654937950870977, 0.14419842183655923, 0.14115839420550458, 0.14047784709222405, 0.14049571462886706, 0.13818879691113792, 0.1379956369831025, 0.13882870569602473, 0.149398324507696, 0.14520869968191066, 0.1313769827201566, 0.14333447289898554, 0.11728672038270234, 0.11659734374312607, 0.11770608713085821, 0.10661103669476844, 0.13294933300463282, 0.1133352805918017, 0.1140511945282083, 0.11058541481581312]
function neid_all_lines_gpu(time_stamps, granulation_status, LD_type, ext_toggle, spot_toggle)

    # NEID location
    obs_lat = 31.9583 
    obs_long = -111.5967  
    alt = 2.097938 

    # set up paramaters for disk
    N = 197
    Nt = length(time_stamps)
    disk = GRASS.Eclipse.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=40, A=14.113, B=-1.698, C=-2.346)

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    depth = GRASS.get_depth(lp)
    lfile = GRASS.get_file(lp)
    #original
    rv = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error = Vector{Vector{Float64}}(undef,size(name)...)
    rv_asymmetric_lsf_inner = Vector{Float64}(undef,size(time_stamps)...)
    rv_error_asymmetric_lsf_inner = Vector{Float64}(undef,size(time_stamps)...)
    resolution = 7e5
    orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
    for i in 1:1#eachindex(lp.λrest) 
        data = jldopen("/storage/home/efg5335/work/completed/NEIDEclipse/data/convolution/NEID_convolution_info$(i).jld2", "r") do file
            Dict(var => read(file, var) for var in variable_names)
        end
        flux_min = data["flux_min"]
        pixel_mean = data["pixel_mean"]
        # set up parameters for synthesis
        lines = [λrest[i]]
        templates = [lfile[i]]
        depths = [depth[i]]
        blueshifts = zeros(length(lines))   # set convective blueshift value

        # make the spec composite type instances
        spec = GRASS.SpecParams(lines=lines, depths=depths,
                        blueshifts=blueshifts, templates=templates, resolution=resolution) 
        
        lambdas, outspec = GRASS.Eclipse.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, obs_long, obs_lat, alt, time_stamps, lines, ext_coeff_array[i], ext_toggle, spot_toggle, true, h5_path; seed_rng=true)
        lambdas .= λ_air_to_vac.(lambdas)
        for j in 1:length(time_stamps)
            outspec_t = outspec[:, j]
            outspec_t = outspec_t ./ maximum(outspec_t)                                                            

            model_pixels, flux_sim, dlam_dpix = neid_asymmetric_lsf.lsf_model_export(orders[i], round(Int, pixel_mean[j]), lambdas, outspec_t, vacwav[i])
            model_wavelength = (dlam_dpix .* model_pixels) .+ vacwav[i]
            v_grid, ccf = GRASS.calc_ccf(model_wavelength, flux_sim, [vacwav[i]], [1.0 - flux_min[j]], 11e4)
            rv_asymmetric_lsf_inner[j], rv_error_asymmetric_lsf_inner[j] = GRASS.calc_rvs_from_ccf(v_grid, ccf)
        end
        rv[i] = rv_asymmetric_lsf_inner
        rv_error[i] = rv_error_asymmetric_lsf_inner
    end
    return rv, rv_error, λrest
end

function neid_october_eclipse_var_on_gpu(LD_type, ext_toggle, spot_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
    rv, rv_error, λrest = deepcopy(neid_all_lines_gpu(neid_october[1:130], true, LD_type, ext_toggle, spot_toggle))

    @save "neid_all_lines_rv_on_$(LD_type)_gpu_disco.jld2" 
    jldopen("neid_all_lines_rv_on_$(LD_type)_gpu_disco.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv_asymmetric_lsf"] = deepcopy(rv) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error)
    end
end

neid_october_eclipse_var_on_gpu("SSD_4parameter", false, true)