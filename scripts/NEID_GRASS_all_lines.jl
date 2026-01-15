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

import PyPlot
plt = PyPlot

using PyCall
py"""
import sys
sys.path.append(".")
"""
neid_symmetric_lsf = pyimport("NEID_LSF")
neid_asymmetric_lsf = pyimport("NeidLsf")
np = pyimport("numpy")

GRASS.get_kernels()

variable_names = ["lambda_min", "lambda_max", "flux_min", "pixel_mean", "pixels", "neid_wavelength"]
# determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)

df_optim_CBOnly = CSV.read("data/Projected_OptimResults_CBOnly.csv", DataFrame; header = false) 
df_optim_CB_MF = CSV.read("data/Projected_OptimResults_CB_MF.csv", DataFrame; header = false)

ext_coeff_array = [0.15452995224327976, 0.15256098077094832, 0.14720055859068512, 0.154895798933504, 0.15181381895180662, 0.15107508233588227, 0.15116772762156633, 0.14882114581650618, 0.14865707189399568, 0.1494903120065096, 0.16011027092744037, 0.15593033972594958, 0.14195968590211427, 0.15401904166429853, 0.1277699772941639, 0.12709315507233226, 0.12820346527304866, 0.11702310600015708, 0.1435320747844216, 0.12380490304619193, 0.12450734135297492, 0.12101777355247835]
function neid_all_lines_gpu(neid_time, granulation_status, LD_type, ext_toggle, model, spot_toggle)
    # convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_time)

    # NEID location
    obs_lat = 31.9583 
    obs_long = -111.5967  
    alt = 2.097938 

    # set up paramaters for disk
    N = 197
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=40)

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    depth = GRASS.get_depth(lp)
    #may want to use optimized depth instead
    lfile = GRASS.get_file(lp)
    #original
    rv = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error = Vector{Vector{Float64}}(undef,size(name)...)
    #NEID symmtetric LSF
    rv_symmetric_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_symmetric_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    #NEID asymmtetric LSF
    rv_asymmetric_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_asymmetric_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    resolution = 7e5
    # loop over lines
    orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
    pixels = [3021:3086, 3062:3128, 2227:2287, 2355:2416, 2464:2526, 2564:2626, 2702:2764, 2738:2801, 2881:2944, 3003:3067, 3057:3093, 3070:3134, 2303:2362, 2529:2589, 3563:3622, 3765:3824, 378:425, 414:461, 478:526, 674:722, 734:781, 802:849]
    for i in eachindex(lp.λrest) 
        if i in vcat(10:22)
            Δv_max=6000.0
        elseif i in vcat(3:9)
            Δv_max=7000.0
        else
            Δv_max=8000.0
        end

        rv_symmetric_lsf_inner = Vector{Float64}(undef,size(time_stamps)...)
        rv_error_symmetric_lsf_inner = Vector{Float64}(undef,size(time_stamps)...)

        rv_asymmetric_lsf_inner = Vector{Float64}(undef,size(time_stamps)...)
        rv_error_asymmetric_lsf_inner = Vector{Float64}(undef,size(time_stamps)...)
        
        println("\t>>> Template: " * string(splitdir(lfile[i])[2]))
        data = jldopen("data/NEID_convolution_info$(i).jld2", "r") do file
            Dict(var => read(file, var) for var in variable_names)
        end
        lambda_min = data["lambda_min"]
        lambda_max = data["lambda_max"]
        flux_min = data["flux_min"]
        pixel_mean = data["pixel_mean"]
        pixels = data["pixels"]
        neid_wavelength = data["neid_wavelength"]

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
            lambdas, outspec = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, obs_long, obs_lat, alt, time_stamps, lines, ext_coeff_array[i], ext_toggle, spot_toggle)
            wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas, outspec, new_res=11e4)
            v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
            rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)
        end
        
        if model == "LD_ext" 
            lambdas, outspec = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, obs_long, obs_lat, alt, time_stamps, lines, ext_coeff_array[i], ext_toggle, spot_toggle)
            wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas, outspec, new_res=11e4)
            v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
            rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)
        end

        if model == "LD_ext_CB" 
            optim_list = df_optim_CBOnly[df_optim_CBOnly.Column1 .== splitext(string(splitdir(lfile[i])[2]))[1], :][1, :]
            lambdas, outspec = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), obs_long, obs_lat, alt, time_stamps, lines, ext_coeff_array[i],
                                    parse(Float64, string(optim_list[2])), parse(Float64, string(optim_list[3])), parse(Float64, string(optim_list[4])), 0.0, 0.0)
            Δv_max=75000.0
            wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas, outspec, new_res=11e4)
            v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4, Δv_max=Δv_max)
            rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)
        end

        if model == "LD_ext_CB_MF" 
            optim_list = df_optim_CB_MF[df_optim_CB_MF.Column1 .== splitext(string(splitdir(lfile[i])[2]))[1], :][1, :]
            lambdas, outspec = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), obs_long, obs_lat, alt, time_stamps, lines, ext_coeff_array[i],
                                    parse(Float64, string(optim_list[2])), parse(Float64, string(optim_list[3])), parse(Float64, string(optim_list[4])), parse(Float64, string(optim_list[5])), parse(Float64, string(optim_list[6])))
            Δv_max=75000.0
            if i == 11
                Δv_max=90000.0
            end
            wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas, outspec, new_res=11e4)
            v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4, Δv_max=Δv_max)
            rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)
        end

        lambdas .= λ_air_to_vac.(lambdas)
        for j in 1:length(time_stamps)
            outspec_t = outspec[:, j]
            outspec_t = outspec_t ./ maximum(outspec_t)                                                            
            # Get the portion of the model for this order
            model_loc_grass = findall(x -> lambda_min[j] ≤ x ≤ lambda_max[j], lambdas)                                                                                              
            wavelength_spacing = np.nanmax(np.diff(lambdas[model_loc_grass]))
            kernel_pixel_arr = np.interp(lambdas[model_loc_grass], neid_wavelength[j], pixels[j])
            kernel_values = pycall(neid_symmetric_lsf.get_variable_lsf_kernel_values, PyObject, pixel_mean[j], kernel_pixel_arr, wavelength_spacing, orders[i]) 
            
            flux_sim = np.convolve(vec(outspec_t), kernel_values, mode="same")
            flux_sim *= wavelength_spacing  
            flux_sim = flux_sim[model_loc_grass] / maximum(flux_sim)  
            wavs_sim = lambdas[model_loc_grass] 
            v_grid, ccf = GRASS.calc_ccf(wavs_sim, flux_sim, [vacwav[i]], [1.0 - flux_min[j]], 11e4, Δv_max=Δv_max)
            rv_symmetric_lsf_inner[j], rv_error_symmetric_lsf_inner[j] = GRASS.calc_rvs_from_ccf(v_grid, ccf)

            model_pixels, flux_sim, dlam_dpix = neid_asymmetric_lsf.lsf_model_export(orders[i], round(Int, pixel_mean[j]), lambdas, outspec_t, vacwav[i])
            model_wavelength = (dlam_dpix .* model_pixels) .+ vacwav[i]
            v_grid, ccf = GRASS.calc_ccf(model_wavelength, flux_sim, [vacwav[i]], [1.0 - flux_min[j]], 11e4, Δv_max=Δv_max)
            rv_asymmetric_lsf_inner[j], rv_error_asymmetric_lsf_inner[j] = GRASS.calc_rvs_from_ccf(v_grid, ccf)
        end

        rv_symmetric_lsf[i] = rv_symmetric_lsf_inner
        rv_error_symmetric_lsf[i] = rv_error_symmetric_lsf_inner

        rv_asymmetric_lsf[i] = rv_asymmetric_lsf_inner
        rv_error_asymmetric_lsf[i] = rv_error_asymmetric_lsf_inner
    end
    return rv, rv_error, λrest, rv_symmetric_lsf, rv_error_symmetric_lsf, rv_asymmetric_lsf, rv_error_asymmetric_lsf
end

function neid_october_eclipse_var_off_gpu(LD_type, ext_toggle, model, spot_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, λrest, rv_symmetric_lsf, rv_error_symmetric_lsf, rv_asymmetric_lsf, rv_error_asymmetric_lsf = deepcopy(neid_all_lines_gpu(neid_october[1:130], false, LD_type, ext_toggle, model, spot_toggle))
    if model == "LD" 
        @save "neid_all_lines_rv_off_$(LD_type)_gpu.jld2"
        jldopen("neid_all_lines_rv_off_$(LD_type)_gpu.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end

    if model == "LD_ext" 
        @save "neid_all_lines_rv_off_$(LD_type)_gpu_ext.jld2" 
        jldopen("neid_all_lines_rv_off_$(LD_type)_gpu_ext.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end

    if model == "LD_ext_CB" 
        @save "neid_all_lines_rv_off_$(LD_type)_gpu_ext_CB_optim.jld2" 
        jldopen("neid_all_lines_rv_off_$(LD_type)_gpu_ext_CB_optim.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end

    if model == "LD_ext_CB_MF" 
        @save "neid_all_lines_rv_off_$(LD_type)_gpu_ext_CB_MF_optim.jld2" 
        jldopen("neid_all_lines_rv_off_$(LD_type)_gpu_ext_CB_MF_optim.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end
end

function neid_october_eclipse_var_on_gpu(LD_type, ext_toggle, model, spot_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, λrest, rv_symmetric_lsf, rv_error_symmetric_lsf, rv_asymmetric_lsf, rv_error_asymmetric_lsf = deepcopy(neid_all_lines_gpu(neid_october[1:130], true, LD_type, ext_toggle, model, spot_toggle))
    if model == "LD" 
        @save "neid_all_lines_rv_on_$(LD_type)_gpu.jld2"
        jldopen("neid_all_lines_rv_on_$(LD_type)_gpu.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end

    if model == "LD_ext" 
        @save "neid_all_lines_rv_on_$(LD_type)_gpu_ext.jld2" 
        jldopen("neid_all_lines_rv_on_$(LD_type)_gpu_ext.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end

    if model == "LD_ext_CB" 
        @save "neid_all_lines_rv_on_$(LD_type)_gpu_ext_CB_optim.jld2" 
        jldopen("neid_all_lines_rv_on_$(LD_type)_gpu_ext_CB_optim.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end

    if model == "LD_ext_CB_MF" 
        @save "neid_all_lines_rv_on_$(LD_type)_gpu_ext_CB_MF_optim.jld2" 
        jldopen("neid_all_lines_rv_on_$(LD_type)_gpu_ext_CB_MF_optim.jld2", "a+") do file
            file["name"] = deepcopy(λrest)
            file["rv"] = deepcopy(rv) 
            file["rv_error"] = deepcopy(rv_error)
            file["rv_symmetric_lsf"] = deepcopy(rv_symmetric_lsf) 
            file["rv_error_symmetric_lsf"] = deepcopy(rv_error_symmetric_lsf)
            file["rv_asymmetric_lsf"] = deepcopy(rv_asymmetric_lsf) 
            file["rv_error_asymmetric_lsf"] = deepcopy(rv_error_asymmetric_lsf)
        end
    end
end

neid_october_eclipse_var_off_gpu("SSD_4parameter", true, "LD_ext_CB_MF", true)