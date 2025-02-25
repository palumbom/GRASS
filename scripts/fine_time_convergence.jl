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
using EchelleCCFs: λ_vac_to_air

using PyCall
py"""
import sys
sys.path.append('.')
"""
neid_lsf = pyimport("NEID_LSF")

GRASS.get_kernels()

function gpu_N_convergence(neid_time, granulation_status, LD_type, ext_toggle)
    # convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_time)

    variable_names = ["wavelength", "variance"]

    # Open the JLD2 file and read the variables into a dictionary
    data = jldopen("data/NEID_October_eclipse_spectrum_info.jld2", "r") do file
        Dict(var => read(file, var) for var in variable_names)
    end

    NEID_wavelength = deepcopy(data["wavelength"])
    NEID_variance = deepcopy(data["variance"])

    # NEID location
    obs_lat = 31.9583 
    obs_long = -111.5967  
    alt = 2.097938 

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    depth = GRASS.get_depth(lp)
    lfile = GRASS.get_file(lp)

    N_range = 50:5:500 #10:5:200
    count = 1
    rv = Vector{Vector{Float64}}(undef,size(N_range)...)
    rv_error = Vector{Vector{Float64}}(undef,size(N_range)...)
    #variance + interpolation
    rv_var = Vector{Vector{Float64}}(undef,size(N_range)...)
    rv_error_var = Vector{Vector{Float64}}(undef,size(N_range)...)
    rvs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    #variance + interpolation + NEID LSF
    rv_var_lsf = Vector{Vector{Float64}}(undef,size(N_range)...)
    rv_error_var_lsf = Vector{Vector{Float64}}(undef,size(N_range)...)
    rvs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    #NEID LSF
    rv_lsf = Vector{Vector{Float64}}(undef,size(N_range)...)
    rv_error_lsf = Vector{Vector{Float64}}(undef,size(N_range)...)
    resolution = 7e5
    # loop over lines
    orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
    pixels = [3021:3086, 3062:3128, 2227:2287, 2355:2416, 2464:2526, 2564:2626, 2702:2764, 2738:2801, 2881:2944, 3003:3067, 3057:3093, 3070:3134, 2303:2362, 2529:2589, 3563:3622, 3765:3824, 378:425, 414:461, 478:526, 674:722, 734:781, 802:849]

    i = 9
    println("\t>>> Template: " * string(splitdir(lfile[i])[2]))
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

    for j in N_range
        # set up paramaters for disk
        N = j
        Nt = length(time_stamps)
        disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=40)

        lambdas_cpu, outspec_cpu = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, 
                                                            obs_long, obs_lat, alt, time_stamps, 
                                                            lines, "three", ext_toggle)
        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv[count], rv_error[count] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var[t], sigs_cpu_var[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
        end

        rv_var[count] = deepcopy(rvs_cpu_var)
        rv_error_var[count] = deepcopy(sigs_cpu_var)

        kernel = pycall(neid_lsf.get_variable_lsf_kernel_values, PyObject, mean(pixels[i]), collect(pixels[i]), maximum(diff(lambdas_cpu)), orders[i])                                                   
        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, kernel, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv_lsf[count], rv_error_lsf[count] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var_lsf[t], sigs_cpu_var_lsf[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
        end

        rv_var_lsf[count] = deepcopy(rvs_cpu_var_lsf)
        rv_error_var_lsf[count] = deepcopy(sigs_cpu_var_lsf)

        count = count + 1
    end
    return rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, N_range
end

function neid_october_eclipse_var_off_gpu_convergence(LD_type, ext_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, N_range = deepcopy(gpu_N_convergence(neid_october, false, LD_type, ext_toggle))
    @save "N_convergence_5434.jld2"
    jldopen("N_convergence_5434.jld2", "a+") do file
        file["N_range"] = deepcopy(N_range) 
        file["rv"] = deepcopy(rv) 
        file["rv_error"] = deepcopy(rv_error)
        file["rv_var"] = deepcopy(rv_var) 
        file["rv_error_var"] = deepcopy(rv_error_var)
        file["rv_lsf"] = deepcopy(rv_lsf) 
        file["rv_error_lsf"] = deepcopy(rv_error_lsf)
        file["rv_var_lsf"] = deepcopy(rv_var_lsf) 
        file["rv_error_var_lsf"] = deepcopy(rv_error_var_lsf)
    end
end

neid_october_eclipse_var_off_gpu_convergence("SSD", false)