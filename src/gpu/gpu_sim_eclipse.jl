function disk_sim_eclipse_gpu(spec::SpecParams{T1}, disk::DiskParams{T1}, 
                               soldata::GPUSolarData{T2}, gpu_allocs::GPUAllocs{T2},
                               eclipse_allocs::EclipseAllocs{T2}, flux_cpu::AA{T1,2},
                               vels::AA{T1,1}, tloop, tloop_init, templates, idx, 
                               obs_long, obs_lat, alt, time_stamps, wavelength;
                               verbose::Bool=false, seed_rng::Bool=false,
                               skip_times::BitVector=falses(disk.Nt)) where {T1<:AF, T2<:AF}

    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)                           
                                
    # parse out GPU allocations
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux

    # parse out rossiter allocations
    lwavgrid = eclipse_allocs.lwavgrid
    rwavgrid = eclipse_allocs.rwavgrid
    allwavs = eclipse_allocs.allwavs
    allints = eclipse_allocs.allints
    bist = eclipse_allocs.bist
    intt = eclipse_allocs.intt
    widt = eclipse_allocs.widt
    ϕc = eclipse_allocs.ϕc
    θc = eclipse_allocs.θc
    μs = eclipse_allocs.μs
    ld = eclipse_allocs.ld
    ext = eclipse_allocs.ext
    dA = eclipse_allocs.dA
    wts = eclipse_allocs.wts
    xyz = eclipse_allocs.xyz
    cbs = eclipse_allocs.cbs
    z_rot = eclipse_allocs.z_rot
    ax_codes = eclipse_allocs.ax_codes
    keys = eclipse_allocs.keys
    dA_total_proj_mean = eclipse_allocs.dA_total_proj_mean
    mean_intensity = eclipse_allocs.mean_intensity
    mean_weight_v_no_cb = eclipse_allocs.mean_weight_v_no_cb
    mean_weight_v_earth_orb = eclipse_allocs.mean_weight_v_earth_orb
    pole_vector_grid = eclipse_allocs.pole_vector_grid
    SP_sun_pos = eclipse_allocs.SP_sun_pos
    SP_sun_vel = eclipse_allocs.SP_sun_vel
    SP_bary = eclipse_allocs.SP_bary
    SP_bary_pos = eclipse_allocs.SP_bary_pos
    SP_bary_vel = eclipse_allocs.SP_bary_vel
    OP_bary = eclipse_allocs.OP_bary
    mu_grid = eclipse_allocs.mu_grid
    projected_velocities_no_cb = eclipse_allocs.projected_velocities_no_cb
    distance = eclipse_allocs.distance
    v_scalar_grid = eclipse_allocs.v_scalar_grid
    v_earth_orb_proj = eclipse_allocs.v_earth_orb_proj

    # alias the input data from GPUSolarData
    disc_mu_gpu = soldata.mu
    disc_ax_gpu = soldata.ax
    lenall_gpu = soldata.len
    cbsall_gpu = soldata.cbs
    bisall_gpu = soldata.bis
    intall_gpu = soldata.int
    widall_gpu = soldata.wid
    depcontrast_gpu = soldata.dep_contrast

    # set number of threads and blocks for len(μ) gpu kernels
    threads1 = 1024
    blocks1 = cld(CUDA.length(μs), prod(threads1))

    # set number of threads and blocks for trimming functions
    threads2 = (4,4,16)
    blocks2 = cld(length(lenall_gpu) * maximum(lenall_gpu) * 100, prod(threads2))

    # set number of threads and blocks for len(μ) * 100 matrix gpu functions
    threads3 = (16,16)
    blocks3 = cld(CUDA.length(μs) * 100, prod(threads3))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads4 = (16,32)
    blocks4 = cld(CUDA.length(μs) * Nλ, prod(threads4))

    threads5 = 1024
    blocks5 = cld(CUDA.length(prof), prod(threads5))

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        bisall_gpu_loop = CUDA.zeros(T2, CUDA.size(bisall_gpu))
        intall_gpu_loop = CUDA.zeros(T2, CUDA.size(intall_gpu))
        widall_gpu_loop = CUDA.zeros(T2, CUDA.size(widall_gpu))
    end

    # # get weighted disk average cbs
    # @cusync sum_wts_og = CUDA.sum(wts)
    # @cusync z_cbs_avg = CUDA.sum(z_cbs .* wts) / sum_wts_og

    # # calculate how much extra shift is needed
    # extra_z = spec.conv_blueshifts .- z_cbs_avg

    # loop over time
    for t in 1:Nt
        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        calc_eclipse_quantities_gpu!(t, time_stamps[t], obs_long, obs_lat, alt, wavelength, disk, gpu_allocs, eclipse_allocs)

#---------------------

        # get sum wts
        @cusync sum_wts = CUDA.sum(wts)

        # # get weighted disk average cbs
        # @cusync sum_wts_og = CUDA.sum(wts)
        # @cusync z_cbs_avg = CUDA.sum(z_cbs .* wts) / sum_wts_og

        # # calculate how much extra shift is needed
        # extra_z = spec.conv_blueshifts .- z_cbs_avg

        # get weighted sum of velocities
        @cusync @inbounds vels[t] = sum(Array(z_rot) .* c_ms .* Array(wts)) ./ sum_wts

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # get a fresh copy of the untrimmed bisector + width data
            @cusync begin
                CUDA.copyto!(bisall_gpu_loop, bisall_gpu)
                CUDA.copyto!(intall_gpu_loop, intall_gpu)
                CUDA.copyto!(widall_gpu_loop, widall_gpu)
            end

            # trim all the bisector data
            @cusync @cuda threads=threads2 blocks=blocks2 trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                             depcontrast_gpu, lenall_gpu,
                                                                             bisall_gpu_loop, intall_gpu_loop,
                                                                             widall_gpu_loop, bisall_gpu,
                                                                             intall_gpu, widall_gpu)

            # assemble line shape on even int grid
            @cusync @cuda threads=threads3 blocks=blocks3 fill_workspaces!(spec.lines[l], spec.variability[l],
                                                                           extra_z[l], tloop, dat_idx,
                                                                           z_rot, z_cbs, lenall_gpu,
                                                                           bisall_gpu_loop, intall_gpu_loop,
                                                                           widall_gpu_loop, allwavs, allints)

            # do the line synthesis, interp back onto wavelength grid
            @cusync @cuda threads=threads4 blocks=blocks4 line_profile_gpu!(prof, μs, wts, λs, allwavs, allints)

            # copy data from GPU to CPU
            @cusync @cuda threads=threads5 blocks=blocks5 apply_line!(t, prof, flux, sum_wts)
        end

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # copy over flux
    @cusync flux_cpu .= Array(flux)

    # make sure nothing is still running on GPU
    CUDA.synchronize()
    return nothing
end
