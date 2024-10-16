using GRASS
using Printf
using Revise
using SPICE
using Statistics
using EchelleCCFs
using BenchmarkTools
using JLD2

GRASS.get_kernels()

case = "Boulder"
# granulation = false

if case == "Boulder"
    boulder_timestamps = ["2023-10-14T13:25:37", "2023-10-14T13:27:01", "2023-10-14T13:28:23", "2023-10-14T13:29:42", "2023-10-14T13:31:08", "2023-10-14T13:32:25", "2023-10-14T13:33:49", "2023-10-14T13:35:10", "2023-10-14T13:36:33", "2023-10-14T13:37:56", "2023-10-14T13:39:18", "2023-10-14T13:40:37", "2023-10-14T13:42:01", "2023-10-14T13:43:22", "2023-10-14T13:44:46", "2023-10-14T13:46:09", "2023-10-14T13:47:26", "2023-10-14T13:48:50", "2023-10-14T13:50:14", "2023-10-14T13:51:38", "2023-10-14T13:52:57", "2023-10-14T13:53:55", "2023-10-14T13:59:00", "2023-10-14T13:59:41", "2023-10-14T14:03:04", "2023-10-14T14:03:20", "2023-10-14T14:05:15", "2023-10-14T14:06:34", "2023-10-14T14:12:22", "2023-10-14T14:13:27", "2023-10-14T14:14:48", "2023-10-14T14:16:12", "2023-10-14T14:17:34", "2023-10-14T14:18:56", "2023-10-14T14:20:20", "2023-10-14T14:21:44", "2023-10-14T14:23:02", "2023-10-14T14:24:24", "2023-10-14T14:25:45", "2023-10-14T14:27:10", "2023-10-14T14:28:32", "2023-10-14T14:29:50", "2023-10-14T14:31:13", "2023-10-14T14:32:34", "2023-10-14T14:33:59", "2023-10-14T14:35:21", "2023-10-14T14:36:18", "2023-10-14T15:38:31", "2023-10-14T15:39:34", "2023-10-14T15:40:55", "2023-10-14T15:42:19", "2023-10-14T15:43:37", "2023-10-14T15:44:42", "2023-10-14T15:46:30", "2023-10-14T15:47:48", "2023-10-14T15:49:10", "2023-10-14T15:50:31", "2023-10-14T15:51:53", "2023-10-14T15:53:14", "2023-10-14T15:54:37", "2023-10-14T15:55:59", "2023-10-14T15:57:21", "2023-10-14T15:58:40", "2023-10-14T16:00:04", "2023-10-14T16:01:28", "2023-10-14T16:02:49", "2023-10-14T16:04:11", "2023-10-14T16:05:32", "2023-10-14T16:06:54", "2023-10-14T16:08:18", "2023-10-14T16:09:39", "2023-10-14T16:11:01", "2023-10-14T16:12:25", "2023-10-14T16:13:46", "2023-10-14T16:15:08", "2023-10-14T16:16:27", "2023-10-14T16:17:50", "2023-10-14T16:19:14", "2023-10-14T16:20:36", "2023-10-14T16:21:57", "2023-10-14T16:23:19", "2023-10-14T16:24:40", "2023-10-14T16:26:04", "2023-10-14T16:27:26", "2023-10-14T16:28:47", "2023-10-14T16:30:09", "2023-10-14T16:31:32", "2023-10-14T16:32:54", "2023-10-14T16:34:16", "2023-10-14T16:35:37", "2023-10-14T16:36:58", "2023-10-14T16:39:52", "2023-10-14T16:41:05", "2023-10-14T16:42:32", "2023-10-14T16:43:51", "2023-10-14T16:45:12", "2023-10-14T16:46:23", "2023-10-14T16:48:03", "2023-10-14T16:49:17", "2023-10-14T16:50:42", "2023-10-14T16:51:58", "2023-10-14T16:53:24", "2023-10-14T16:54:45", "2023-10-14T16:56:06", "2023-10-14T16:57:31", "2023-10-14T16:58:53", "2023-10-14T17:00:13", "2023-10-14T17:01:35", "2023-10-14T17:02:56", "2023-10-14T17:04:20", "2023-10-14T17:05:42", "2023-10-14T17:07:03", "2023-10-14T17:08:27", "2023-10-14T17:09:49", "2023-10-14T17:11:10", "2023-10-14T17:12:32", "2023-10-14T17:13:53", "2023-10-14T17:15:17", "2023-10-14T17:16:39", "2023-10-14T17:18:00", "2023-10-14T17:19:21", "2023-10-14T17:20:43", "2023-10-14T17:22:07", "2023-10-14T17:23:28", "2023-10-14T17:24:50", "2023-10-14T17:26:12", "2023-10-14T17:27:33", "2023-10-14T17:28:41", "2023-10-14T18:41:52", "2023-10-14T18:42:46", "2023-10-14T18:44:08", "2023-10-14T18:45:29", "2023-10-14T18:46:50", "2023-10-14T18:48:12", "2023-10-14T18:49:36", "2023-10-14T18:50:39", "2023-10-14T18:55:41", "2023-10-14T18:56:25", "2023-10-14T18:57:51", "2023-10-14T18:59:09", "2023-10-14T19:00:33", "2023-10-14T19:01:54", "2023-10-14T19:03:13", "2023-10-14T19:04:37", "2023-10-14T20:37:30", "2023-10-14T23:07:01", "2023-10-14T23:07:58", "2023-10-14T23:09:18"]
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(boulder_timestamps)

    variable_names = ["zenith", "dA_total_proj", "idx1", "idx3", "mu_grid", "z_rot_sub", "mu", "ax_codes", "dA", "N", "dA_total_proj_mean", "mean_weight_v_no_cb", "mean_weight_v_earth_orb"]
    # Open the JLD2 file and read the variables into a dictionary
    data = jldopen("data/solar_disk/boulder_october_N_50.jld2", "r") do file
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
    dA_total_proj_mean = deepcopy(data["dA_total_proj_mean"])
    mean_weight_v_no_cb = deepcopy(data["mean_weight_v_no_cb"])
    mean_weight_v_earth_orb = deepcopy(data["mean_weight_v_earth_orb"])

    #Boulder location
    obs_lat = 39.995380
    obs_long = -105.262390
    alt = 1.6523

    # set up paramaters for disk
    N = 50
    Nt = length(time_stamps)

    # set up parameters for spectrum
    lines = [15652.79] # array of line centers 
    # depths = [0.25]   # array of line depths
    # templates = ["FeI_5434"] # template data to use
    # if granulation == true
    #     variability = trues(length(lines))  # whether or not the bisectors should "dance"
    # end
    # if granulation == false
    #     variability = falses(length(lines))
    # end
    # blueshifts = zeros(length(lines))   # set convective blueshift value
    # resolution = 7e5                    # spectral resolution
    # new_res = 1e6

    # # make the disk and spec composite type instances
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)
    wsp = GRASS.SynthWorkspaceEclipse(disk, 1, Nt, verbose=true)
    # spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
    #                 blueshifts=blueshifts, templates=templates, resolution=resolution)  
    
    LD_arr = ["HD", "K300", "KSSD", "NIR"]
    # rv = Vector{Vector{Float64}}(undef,size(LD_arr)...)
    # rv_error = Vector{Vector{Float64}}(undef,size(LD_arr)...)
    RV_list_no_cb_final = Vector{Vector{Float64}}(undef,size(LD_arr)...)
    intensity_list_final = Vector{Vector{Float64}}(undef,size(LD_arr)...)
    for i in 1:length(LD_arr)
        # # actually synthesize the spectra
        # lambdas_cpu, outspec_cpu = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_arr[i], zenith_mean,
        #                                     dA_total_proj, idx1, idx3, mu_grid, z_rot_sub,
        #                                     stored_μs, stored_ax_codes, stored_dA, "three", ext_toggle = false, verbose=true, use_gpu=false)
        # wavs_sim, flux_sim = GRASS.convolve_gauss(lambdas_cpu, outspec_cpu, new_res=new_res, oversampling=4.0)

        # #measure velocities
        # v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wavs_sim, flux_sim, spec)
        # rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        # rv[i] = rvs_cpu
        # rv_error[i] = sigs_cpu
        
        RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
        intensity_list = Vector{Float64}(undef,size(time_stamps)...)
        for t in 1:length(time_stamps)
            mean_intensity = GRASS.eclipse_compute_intensity(disk, lines, [2.0], LD_arr[i], idx1[t], idx3[t],
                                mu_grid[t], mean_weight_v_no_cb[:, :, t], mean_weight_v_earth_orb[:, :, t],
                                z_rot_sub[t], dA_total_proj[t], wsp.ld, wsp.z_rot, zenith_mean[t], false, wsp.ext)
            idx_grid = mean_intensity[:, :, 1] .> 0.0

            brightness = mean_intensity[:, :, 1] .* dA_total_proj_mean[:, :, t]

            cheapflux = sum(view(brightness, idx_grid))

            #determine final mean weighted velocity for disk grid
            final_weight_v_no_cb = sum(view(mean_weight_v_no_cb[:, :, t] .* brightness, idx_grid)) / cheapflux 
            final_weight_v_no_cb += mean(view(mean_weight_v_earth_orb[:, :, t], idx_grid))  

            RV_list_no_cb[t] = deepcopy(final_weight_v_no_cb)
            intensity_list[t] = deepcopy(cheapflux)
        end

        RV_list_no_cb_final[i] = deepcopy(RV_list_no_cb)
        intensity_list_final[i] = deepcopy(intensity_list)
    end

    # @save "boulder_rv_off.jld2"
    # jldopen("boulder_rv_off.jld2", "a+") do file
    #     file["rv"] = deepcopy(rv) 
    #     file["rv_error"] = deepcopy(rv_error)
    # end

    @save "boulder_rv_projected.jld2"
        jldopen("boulder_rv_projected.jld2", "a+") do file
        file["RV_list_no_cb"] = deepcopy(RV_list_no_cb_final) 
        file["intensity_list"] = deepcopy(intensity_list_final)
    end
end

########################################################################
#TO BE CLEAR
#grass = src: gpu (finalize and clear) + structures (GPU ones need to be finalized)
#grass once gpu done run with gpu + clear loose in src and loose outside src - that follow GPU pipeline

# if case == "GPU_test" 
#     epoch = utc2et.("2023-10-14T15:26:45")

#     #NEID location
#     obs_lat = 31.9583 
#     obs_long = -111.5967  
#     alt = 2.097938 

#     # set up paramaters for disk
#     N = 2
#     Nsubgrid=2
#     Nt = length(epoch)
    
#     # set up parameters for spectrum
#     lines = [6173.0, 6173.0] # array of line centers 
#     depths = [0.6, 0.6]   # array of line depths
#     templates = ["FeI_6173", "FeI_6173"] # template data to use
        # if granulation == true
        #     variability = trues(length(lines))  # whether or not the bisectors should "dance"
        # end
        # if granulation == false
        #     variability = falses(length(lines))
        # end
#     blueshifts = zeros(length(lines))   # set convective blueshift value
#     resolution = 7e5                    # spectral resolution

#     # make the disk and spec composite type instances
#     disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=Nsubgrid)
#     spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
#                     blueshifts=blueshifts, templates=templates, resolution=resolution) 
    
#     gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, Int(length(lines)))
#     GRASS.calc_eclipse_quantities_gpu!(epoch, obs_long, obs_lat, alt, lines./10.0, disk, gpu_allocs)

#     println(Array(gpu_allocs.ld)) #plt.imshow(diff between two arrays)
#     #maximum(abs.(wsp.ld[:,:,1] .- Array(gpu_allocs.ld)[:,:,1]))
#     #plt.imshow(wsp.ld[:,:,1] .- Array(gpu_allocs.ld)[:,:,1])
#     println("----------------------")
#     wsp = GRASS.SynthWorkspaceEclipse(disk, Int(length(lines)))
#     mem = GRASS.GeoWorkspaceEclipse(disk, Int(length(lines)))
#     GRASS.eclipse_compute_quantities!(disk, epoch, obs_long, obs_lat, alt, lines./10.0, wsp.ϕc, wsp.θc, 
#                                                 wsp.μs, wsp.ld, wsp.ext, wsp.dA, wsp.xyz, wsp.wts, wsp.z_rot, wsp.ax_codes,
#                                                 mem.dA_total_proj_mean, mem.mean_intensity, mem.mean_weight_v_no_cb,
#                                                 mem.mean_weight_v_earth_orb, mem.pole_vector_grid,
#                                                 mem.SP_sun_pos, mem.SP_sun_vel, mem.SP_bary, mem.SP_bary_pos,
#                                                 mem.SP_bary_vel, mem.OP_bary, mem.mu_grid, mem.projected_velocities_no_cb, 
#                                                 mem.distance, mem.v_scalar_grid, mem.v_earth_orb_proj)

#     println(wsp.ld)       
# end