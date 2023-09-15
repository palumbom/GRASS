function disk_sim_gpu(spec::SpecParams{T1}, disk::DiskParams{T1}, soldata::GPUSolarData{T2},
                      gpu_allocs::GPUAllocs{T2}, flux::AA{T1,2}; verbose::Bool=false,
                      seed_rng::Bool=false,  skip_times::BitVector=falses(disk.Nt)) where {T1<:AF, T2<:AF}
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # parse out composite type
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof

    μs = gpu_allocs.μs
    wts = gpu_allocs.wts
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs

    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints

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

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        bisall_gpu_loop = CUDA.zeros(T2, CUDA.size(bisall_gpu))
        intall_gpu_loop = CUDA.zeros(T2, CUDA.size(intall_gpu))
        widall_gpu_loop = CUDA.zeros(T2, CUDA.size(widall_gpu))
    end

    # get weighted disk average cbs
    @cusync sum_wts = CUDA.sum(wts)
    @cusync z_cbs_avg = CUDA.sum(z_cbs .* wts) / sum_wts

    # calculate how much extra shift is needed
    extra_z = spec.conv_blueshifts .- z_cbs_avg

    # loop over time
    for t in 1:Nt
        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # re-zero the line profile holder
            @cusync prof .= zero(T2)

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
            @cusync @inbounds flux[:,t] .*= Array(prof) ./ sum_wts
        end

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # make sure nothing is still running on GPU
    CUDA.synchronize()
    return nothing
end

function disk_sim_rossiter_gpu(spec::SpecParams{T1}, disk::DiskParams{T1}, planet::Planet{T1},
                               soldata::GPUSolarData{T2}, gpu_allocs::GPUAllocs{T2},
                               ros_allocs::RossiterAllocsGPU{T2}, flux::AA{T1,2},
                               vels::AA{T1,1}; verbose::Bool=false, seed_rng::Bool=false,
                               skip_times::BitVector=falses(disk.Nt)) where {T1<:AF, T2<:AF}
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # parse out GPU allocations
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof

    z_cbs = gpu_allocs.z_cbs

    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints

    # parse out rossiter allocations
    μs = ros_allocs.μs
    wts = ros_allocs.wts
    z_rot = ros_allocs.z_rot
    xyz_planet = ros_allocs.xyz_planet
    xyz_star = ros_allocs.xyz_star

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

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        bisall_gpu_loop = CUDA.zeros(T2, CUDA.size(bisall_gpu))
        intall_gpu_loop = CUDA.zeros(T2, CUDA.size(intall_gpu))
        widall_gpu_loop = CUDA.zeros(T2, CUDA.size(widall_gpu))
    end

    # get weighted disk average cbs
    @cusync sum_wts_og = CUDA.sum(wts)
    @cusync z_cbs_avg = CUDA.sum(z_cbs .* wts) / sum_wts_og

    # calculate how much extra shift is needed
    extra_z = spec.conv_blueshifts .- z_cbs_avg

    # loop over time
    for t in 1:Nt
        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        # recopy weights, etc. from unobstructed disk
        @cusync begin
            CUDA.copyto!(ros_allocs.μs, gpu_allocs.μs)
            CUDA.copyto!(ros_allocs.wts, gpu_allocs.wts)
            CUDA.copyto!(ros_allocs.z_rot, gpu_allocs.z_rot)
        end

        # TODO: this isn't needed if planet is out of transit at time step
        # re-calculate patch weights, etc. for occulted patches
        calc_rossiter_quantities_gpu!(xyz_planet, t, planet, disk, gpu_allocs, ros_allocs)

        # get sum wts
        @cusync sum_wts = CUDA.sum(wts)

        # get weighted sum of velocities
        @cusync @inbounds vels[t] = sum(Array(z_rot) .* c_ms .* Array(wts)) ./ sum_wts

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # re-zero the line profile holder
            @cusync prof .= zero(T2)

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
            @cusync @inbounds flux[:,t] .*= Array(prof) ./ sum_wts
        end

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # make sure nothing is still running on GPU
    CUDA.synchronize()
    return nothing
end
