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

function neid_all_lines(neid_time, granulation_status, filename, LD_type, ext_toggle)
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

    variable_names = ["zenith", "dA_total_proj", "idx1", "idx3", "mu_grid", "z_rot_sub", "mu", "ax_codes", "dA", "N"]

    # Open the JLD2 file and read the variables into a dictionary
    data = jldopen("data/solar_disk/$(filename).jld2", "r") do file
        Dict(var => read(file, var) for var in variable_names)
    end

    zenith_mean = deepcopy(data["zenith"])
    dA_total_proj = deepcopy(data["dA_total_proj"])
    idx1 = deepcopy(data["idx1"])
    idx3 = deepcopy(data["idx3"])
    mu_grid = deepcopy(data["mu_grid"])
    z_rot_sub = deepcopy(data["z_rot_sub"])
    stored_μs = deepcopy(data["mu"])
    stored_ax_codes = deepcopy(data["ax_codes"])
    stored_dA = deepcopy(data["dA"])

    # set up paramaters for disk
    N = data["N"]
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    depth = GRASS.get_depth(lp)
    lfile = GRASS.get_file(lp)

    #original
    rv = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error = Vector{Vector{Float64}}(undef,size(name)...)
    #variance + interpolation
    rv_var = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_var = Vector{Vector{Float64}}(undef,size(name)...)
    rvs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    #variance + interpolation + NEID LSF
    rv_var_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_var_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rvs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    #NEID LSF
    rv_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    resolution = 7e5
    # loop over lines
    orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
    pixels = [3021:3086, 3062:3128, 2227:2287, 2355:2416, 2464:2526, 2564:2626, 2702:2764, 2738:2801, 2881:2944, 3003:3067, 3057:3093, 3070:3134, 2303:2362, 2529:2589, 3563:3622, 3765:3824, 378:425, 414:461, 478:526, 674:722, 734:781, 802:849]
    for i in eachindex(lp.λrest) 
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

        lambdas_cpu, outspec_cpu = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_type, obs_long, obs_lat, alt, 
                                                            time_stamps, zenith_mean, dA_total_proj, idx1, idx3, 
                                                            mu_grid, z_rot_sub, stored_μs, stored_ax_codes, stored_dA, "three", 
                                                            ext_toggle = ext_toggle, verbose=true, use_gpu=false)
        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var[t], sigs_cpu_var[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
        end

        rv_var[i] = deepcopy(rvs_cpu_var)
        rv_error_var[i] = deepcopy(sigs_cpu_var)

        kernel = pycall(neid_lsf.get_variable_lsf_kernel_values, PyObject, mean(pixels[i]), collect(pixels[i]), maximum(diff(lambdas_cpu)), orders[i])                                                   
        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, kernel, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv_lsf[i], rv_error_lsf[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var_lsf[t], sigs_cpu_var_lsf[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
        end

        rv_var_lsf[i] = deepcopy(rvs_cpu_var_lsf)
        rv_error_var_lsf[i] = deepcopy(sigs_cpu_var_lsf)
    end
    return rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, λrest
end

function neid_all_lines_gpu(neid_time, granulation_status, filename, LD_type, ext_toggle)
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

    # set up paramaters for disk
    N = 50
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    depth = GRASS.get_depth(lp)
    lfile = GRASS.get_file(lp)

    #original
    rv = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error = Vector{Vector{Float64}}(undef,size(name)...)
    #variance + interpolation
    rv_var = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_var = Vector{Vector{Float64}}(undef,size(name)...)
    rvs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    #variance + interpolation + NEID LSF
    rv_var_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_var_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rvs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    #NEID LSF
    rv_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    resolution = 7e5
    # loop over lines
    orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
    pixels = [3021:3086, 3062:3128, 2227:2287, 2355:2416, 2464:2526, 2564:2626, 2702:2764, 2738:2801, 2881:2944, 3003:3067, 3057:3093, 3070:3134, 2303:2362, 2529:2589, 3563:3622, 3765:3824, 378:425, 414:461, 478:526, 674:722, 734:781, 802:849]
    for i in eachindex(lp.λrest) 
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

        lambdas_cpu, unbinned_flux, brightness = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, 
                                                            obs_long, obs_lat, alt, time_stamps, 
                                                            lines, "three", ext_toggle) 
  
        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var[t], sigs_cpu_var[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
        end
        rv_var[i] = deepcopy(rvs_cpu_var)
        rv_error_var[i] = deepcopy(sigs_cpu_var)

        kernel = pycall(neid_lsf.get_variable_lsf_kernel_values, PyObject, mean(pixels[i]), collect(pixels[i]), maximum(diff(lambdas_cpu)), orders[i])                                                   
        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, kernel, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv_lsf[i], rv_error_lsf[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var_lsf[t], sigs_cpu_var_lsf[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
        end

        rv_var_lsf[i] = deepcopy(rvs_cpu_var_lsf)
        rv_error_var_lsf[i] = deepcopy(sigs_cpu_var_lsf)
    end
    return rv, rv_error, λrest, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf
end

function neid_all_lines_gpu_fine_time(neid_time, granulation_status, filename, LD_type, ext_toggle)
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

    # set up paramaters for disk
    N = 50
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    depth = GRASS.get_depth(lp)
    lfile = GRASS.get_file(lp)

    #original
    # brightness = Vector{Vector{Float64}}(undef,size([λrest[1]])...)
    rv = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error = Vector{Vector{Float64}}(undef,size(name)...)
    #variance + interpolation
    rv_var = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_var = Vector{Vector{Float64}}(undef,size(name)...)
    rvs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var = Vector{Float64}(undef,size(time_stamps)...)
    #variance + interpolation + NEID LSF
    rv_var_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_var_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rvs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    sigs_cpu_var_lsf = Vector{Float64}(undef,size(time_stamps)...)
    #NEID LSF
    rv_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    rv_error_lsf = Vector{Vector{Float64}}(undef,size(name)...)
    resolution = 7e5
    # loop over lines
    orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
    pixels = [3021:3086, 3062:3128, 2227:2287, 2355:2416, 2464:2526, 2564:2626, 2702:2764, 2738:2801, 2881:2944, 3003:3067, 3057:3093, 3070:3134, 2303:2362, 2529:2589, 3563:3622, 3765:3824, 378:425, 414:461, 478:526, 674:722, 734:781, 802:849]
    for i in eachindex(lp.λrest) 
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

        lambdas_cpu, unbinned_flux, brightness = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, 
                                                            obs_long, obs_lat, alt, time_stamps, 
                                                            lines, "three", ext_toggle) 
        block_size = 6#56
        # Number of blocks
        num_blocks = size(unbinned_flux, 2) ÷ block_size
        # Initialize a new matrix to store the means
        outspec_cpu = zeros(size(unbinned_flux, 1), num_blocks)
        # Loop through each block and compute the mean of the block
        for i in 1:num_blocks
            # Slice the block (5 columns at a time)
            start_col = (i - 1) * block_size + 1
            end_col = i * block_size
            # outspec_cpu[:, i] = mean(unbinned_flux[:, start_col:end_col], dims=2) 
            outspec_cpu[:, i] = mean(unbinned_flux[:, start_col:end_col][:, 1:4], dims=2) 
        end      

        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv[i], rv_error[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        # NEID_ind = 1
        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var[t], sigs_cpu_var[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
            # if t == (NEID_ind * 56)
            #     NEID_ind += 1
            # end
        end
        rv_var[i] = deepcopy(rvs_cpu_var)
        rv_error_var[i] = deepcopy(sigs_cpu_var)

        kernel = pycall(neid_lsf.get_variable_lsf_kernel_values, PyObject, mean(pixels[i]), collect(pixels[i]), maximum(diff(lambdas_cpu)), orders[i])                                                   
        wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, kernel, new_res=11e4)
        v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths, 11e4)
        rv_lsf[i], rv_error_lsf[i] = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        for t in 1:size(flux_sim, 2)
            itp = interpolate((wavs_sim,), flux_sim[:,t], Gridded(Linear()))
            wavs_sim_itp = λ_vac_to_air.(NEID_wavelength[i][t])
            flux_sim_itp = itp.(wavs_sim_itp)
            # measure velocities
            v_grid_cpu, ccf_cpu, ccf_var_out = (GRASS.calc_ccf(wavs_sim_itp, flux_sim_itp, filter(!=(0), NEID_variance[i][:,t]), lines, depths, 11e4))
            rvs_cpu_var_lsf[t], sigs_cpu_var_lsf[t] = (GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out))
        end

        rv_var_lsf[i] = deepcopy(rvs_cpu_var_lsf)
        rv_error_var_lsf[i] = deepcopy(sigs_cpu_var_lsf)
    end
    return rv, rv_error, λrest, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf
end

function gpu_N_convergence(neid_time, granulation_status, filename, LD_type, ext_toggle)
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

    N_range = 10:5:200 #50:5:500
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

    i = 1
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
        N = 197
        Nt = length(time_stamps)
        disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=j)

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

function neid_october_eclipse_var_off(LD_type, ext_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, λrest = deepcopy(neid_all_lines(neid_october, false, "neid_october_N_50", LD_type, ext_toggle))
    @save "neid_all_lines_rv_off_$(LD_type).jld2"
    jldopen("neid_all_lines_rv_off_$(LD_type).jld2", "a+") do file
        file["name"] = deepcopy(λrest) 
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

function neid_october_eclipse_var_on(LD_type, ext_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, λrest = deepcopy(neid_all_lines(neid_october, true, "neid_october_N_50", LD_type, ext_toggle))
    @save "neid_all_lines_rv_regular_$(LD_type).jld2"
    jldopen("neid_all_lines_rv_regular_$(LD_type).jld2", "a+") do file
        file["name"] = deepcopy(λrest) 
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

function neid_october_eclipse_var_off_gpu(LD_type, ext_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, λrest, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf = deepcopy(neid_all_lines_gpu(neid_october, false, "neid_october_N_50", LD_type, ext_toggle))
    @save "neid_all_lines_rv_off_$(LD_type)_gpu.jld2"
    jldopen("neid_all_lines_rv_off_$(LD_type)_gpu.jld2", "a+") do file
        file["name"] = deepcopy(λrest)
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

function neid_october_eclipse_var_on_gpu(LD_type, ext_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, λrest = deepcopy(neid_all_lines_gpu(neid_october, true, "neid_october_N_50", LD_type, ext_toggle))
    @save "neid_all_lines_rv_regular_$(LD_type)_gpu.jld2"
    jldopen("neid_all_lines_rv_regular_$(LD_type)_gpu.jld2", "a+") do file
        file["name"] = deepcopy(λrest) 
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

function neid_october_nxt_day_var_off(LD_type, ext_toggle)
    neid_1015 = ["2023-10-15T16:33:02.500", "2023-10-15T16:34:25.500", "2023-10-15T16:35:48.500", "2023-10-15T16:37:10.500", "2023-10-15T16:38:33.500", "2023-10-15T16:39:56.500", "2023-10-15T16:41:18.500", "2023-10-15T16:42:41.500", "2023-10-15T16:44:04.500", "2023-10-15T16:45:26.500", "2023-10-15T16:46:49.500", "2023-10-15T16:48:12.500", "2023-10-15T16:49:34.500", "2023-10-15T16:50:57.500", "2023-10-15T16:52:20.500", "2023-10-15T16:53:42.500", "2023-10-15T16:55:05.500", "2023-10-15T16:56:28.500", "2023-10-15T16:57:50.500", "2023-10-15T16:59:13.500", "2023-10-15T17:00:36.500", "2023-10-15T17:01:59.500", "2023-10-15T17:03:21.500", "2023-10-15T17:04:44.500", "2023-10-15T17:06:07.500", "2023-10-15T17:07:29.500", "2023-10-15T17:08:52.500", "2023-10-15T17:10:15.500", "2023-10-15T17:11:37.500", "2023-10-15T17:13:00.500", "2023-10-15T17:14:23.500", "2023-10-15T17:15:45.500", "2023-10-15T17:17:08.500", "2023-10-15T17:18:31.500", "2023-10-15T17:19:53.500", "2023-10-15T17:21:16.500", "2023-10-15T17:22:39.500", "2023-10-15T17:24:01.500", "2023-10-15T17:25:24.500", "2023-10-15T17:26:47.500", "2023-10-15T17:28:09.500", "2023-10-15T17:29:32.500", "2023-10-15T17:30:55.500", "2023-10-15T17:32:18.500", "2023-10-15T17:33:40.500", "2023-10-15T17:35:03.500", "2023-10-15T17:36:26.500", "2023-10-15T17:37:48.500", "2023-10-15T17:39:11.500", "2023-10-15T17:40:34.500", "2023-10-15T17:41:56.500", "2023-10-15T17:43:19.500", "2023-10-15T17:44:42.500", "2023-10-15T17:46:04.500", "2023-10-15T17:47:27.500", "2023-10-15T17:48:50.500", "2023-10-15T17:50:12.500", "2023-10-15T17:51:35.500", "2023-10-15T17:52:58.500", "2023-10-15T17:54:20.500", "2023-10-15T17:57:10.500", "2023-10-15T17:58:33.500", "2023-10-15T17:59:55.500", "2023-10-15T18:01:18.500", "2023-10-15T18:02:41.500", "2023-10-15T18:04:03.500", "2023-10-15T18:05:26.500", "2023-10-15T18:06:49.500", "2023-10-15T18:08:11.500", "2023-10-15T18:09:34.500", "2023-10-15T18:10:57.500", "2023-10-15T18:12:19.500", "2023-10-15T18:13:42.500", "2023-10-15T18:15:05.500", "2023-10-15T18:16:27.500", "2023-10-15T18:17:50.500", "2023-10-15T18:19:13.500", "2023-10-15T18:20:35.500", "2023-10-15T18:21:58.500", "2023-10-15T18:23:21.500", "2023-10-15T18:24:44.500", "2023-10-15T18:26:06.500", "2023-10-15T18:27:29.500", "2023-10-15T18:28:52.500", "2023-10-15T18:30:14.500", "2023-10-15T18:31:37.500", "2023-10-15T18:33:00.500", "2023-10-15T18:34:22.500", "2023-10-15T18:35:45.500", "2023-10-15T18:37:08.500", "2023-10-15T18:38:30.500", "2023-10-15T18:39:53.500", "2023-10-15T18:41:16.500", "2023-10-15T18:42:38.500", "2023-10-15T18:44:01.500", "2023-10-15T18:45:24.500", "2023-10-15T18:46:46.500", "2023-10-15T18:48:09.500", "2023-10-15T18:49:32.500", "2023-10-15T18:50:54.500", "2023-10-15T18:52:17.500", "2023-10-15T18:53:40.500", "2023-10-15T18:55:03.500", "2023-10-15T18:56:25.500", "2023-10-15T18:57:48.500", "2023-10-15T18:59:11.500", "2023-10-15T19:00:33.500", "2023-10-15T19:01:56.500", "2023-10-15T19:03:19.500", "2023-10-15T19:04:41.500", "2023-10-15T19:06:04.500"]

    rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, λrest = deepcopy(neid_all_lines(neid_1015, false, "neid_october_N_50_nxt_day", LD_type, ext_toggle))
    @save "neid_all_lines_rv_off_$(LD_type)_nxt_day.jld2"
    jldopen("neid_all_lines_rv_off_$(LD_type)_nxt_day.jld2", "a+") do file
        file["name"] = deepcopy(λrest) 
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

function neid_october_nxt_day_var_on(LD_type, ext_toggle)
    neid_1015 = ["2023-10-15T16:33:02.500", "2023-10-15T16:34:25.500", "2023-10-15T16:35:48.500", "2023-10-15T16:37:10.500", "2023-10-15T16:38:33.500", "2023-10-15T16:39:56.500", "2023-10-15T16:41:18.500", "2023-10-15T16:42:41.500", "2023-10-15T16:44:04.500", "2023-10-15T16:45:26.500", "2023-10-15T16:46:49.500", "2023-10-15T16:48:12.500", "2023-10-15T16:49:34.500", "2023-10-15T16:50:57.500", "2023-10-15T16:52:20.500", "2023-10-15T16:53:42.500", "2023-10-15T16:55:05.500", "2023-10-15T16:56:28.500", "2023-10-15T16:57:50.500", "2023-10-15T16:59:13.500", "2023-10-15T17:00:36.500", "2023-10-15T17:01:59.500", "2023-10-15T17:03:21.500", "2023-10-15T17:04:44.500", "2023-10-15T17:06:07.500", "2023-10-15T17:07:29.500", "2023-10-15T17:08:52.500", "2023-10-15T17:10:15.500", "2023-10-15T17:11:37.500", "2023-10-15T17:13:00.500", "2023-10-15T17:14:23.500", "2023-10-15T17:15:45.500", "2023-10-15T17:17:08.500", "2023-10-15T17:18:31.500", "2023-10-15T17:19:53.500", "2023-10-15T17:21:16.500", "2023-10-15T17:22:39.500", "2023-10-15T17:24:01.500", "2023-10-15T17:25:24.500", "2023-10-15T17:26:47.500", "2023-10-15T17:28:09.500", "2023-10-15T17:29:32.500", "2023-10-15T17:30:55.500", "2023-10-15T17:32:18.500", "2023-10-15T17:33:40.500", "2023-10-15T17:35:03.500", "2023-10-15T17:36:26.500", "2023-10-15T17:37:48.500", "2023-10-15T17:39:11.500", "2023-10-15T17:40:34.500", "2023-10-15T17:41:56.500", "2023-10-15T17:43:19.500", "2023-10-15T17:44:42.500", "2023-10-15T17:46:04.500", "2023-10-15T17:47:27.500", "2023-10-15T17:48:50.500", "2023-10-15T17:50:12.500", "2023-10-15T17:51:35.500", "2023-10-15T17:52:58.500", "2023-10-15T17:54:20.500", "2023-10-15T17:57:10.500", "2023-10-15T17:58:33.500", "2023-10-15T17:59:55.500", "2023-10-15T18:01:18.500", "2023-10-15T18:02:41.500", "2023-10-15T18:04:03.500", "2023-10-15T18:05:26.500", "2023-10-15T18:06:49.500", "2023-10-15T18:08:11.500", "2023-10-15T18:09:34.500", "2023-10-15T18:10:57.500", "2023-10-15T18:12:19.500", "2023-10-15T18:13:42.500", "2023-10-15T18:15:05.500", "2023-10-15T18:16:27.500", "2023-10-15T18:17:50.500", "2023-10-15T18:19:13.500", "2023-10-15T18:20:35.500", "2023-10-15T18:21:58.500", "2023-10-15T18:23:21.500", "2023-10-15T18:24:44.500", "2023-10-15T18:26:06.500", "2023-10-15T18:27:29.500", "2023-10-15T18:28:52.500", "2023-10-15T18:30:14.500", "2023-10-15T18:31:37.500", "2023-10-15T18:33:00.500", "2023-10-15T18:34:22.500", "2023-10-15T18:35:45.500", "2023-10-15T18:37:08.500", "2023-10-15T18:38:30.500", "2023-10-15T18:39:53.500", "2023-10-15T18:41:16.500", "2023-10-15T18:42:38.500", "2023-10-15T18:44:01.500", "2023-10-15T18:45:24.500", "2023-10-15T18:46:46.500", "2023-10-15T18:48:09.500", "2023-10-15T18:49:32.500", "2023-10-15T18:50:54.500", "2023-10-15T18:52:17.500", "2023-10-15T18:53:40.500", "2023-10-15T18:55:03.500", "2023-10-15T18:56:25.500", "2023-10-15T18:57:48.500", "2023-10-15T18:59:11.500", "2023-10-15T19:00:33.500", "2023-10-15T19:01:56.500", "2023-10-15T19:03:19.500", "2023-10-15T19:04:41.500", "2023-10-15T19:06:04.500"]

    rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, λrest = deepcopy(neid_all_lines(neid_1015, true, "neid_october_N_50_nxt_day", LD_type, ext_toggle))
    @save "neid_all_lines_rv_regular_$(LD_type)_nxt_day.jld2"
    jldopen("neid_all_lines_rv_regular_$(LD_type)_nxt_day.jld2", "a+") do file
        file["name"] = deepcopy(λrest) 
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

function neid_october_eclipse_var_off_gpu_fine_time(LD_type, ext_toggle)
    df = CSV.File("data/NEID_October_15s_time.csv") |> DataFrame

    # Display the DataFrame
    neid_october = string.(df[!, "FineTime"])

    rv, rv_error, λrest, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf = deepcopy(neid_all_lines_gpu_fine_time(neid_october, false, "neid_october_N_50", LD_type, ext_toggle))
    @save "neid_all_lines_rv_off_$(LD_type)_gpu_15s_time.jld2"
    jldopen("neid_all_lines_rv_off_$(LD_type)_gpu_15s_time.jld2", "a+") do file
        file["name"] = deepcopy(λrest)
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

function neid_october_eclipse_var_off_gpu_convergence(LD_type, ext_toggle)
    neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]

    rv, rv_error, rv_var, rv_error_var, rv_var_lsf, rv_error_var_lsf, rv_lsf, rv_error_lsf, N_range = deepcopy(gpu_N_convergence(neid_october, false, "neid_october_N_50", LD_type, ext_toggle))
    @save "subgrid_convergence.jld2"
    jldopen("subgrid_convergence.jld2", "a+") do file
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

neid_october_eclipse_var_off_gpu_fine_time("SSD", true)