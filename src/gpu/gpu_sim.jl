function disk_sim_gpu(spec::SpecParams{T1}, disk::DiskParams{T1}, soldata::GPUSolarData{T2},
                      gpu_allocs::GPUAllocs{T2}, flux_cpu::AA{T1,3}; verbose::Bool=false,
                      seed_rng::Bool=false,  skip_times::BitVector=falses(disk.Nt),
                      show_progress::Bool=true) where {T1<:AF, T2<:AF}
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # parse out composite type
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux

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

    # allocate destinations for interpolations
    @cusync begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
    end

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

    # get weighted disk average cbs
    @cusync sum_wts = CUDA.sum(wts)
    @cusync z_cbs_avg = CUDA.sum(z_cbs .* wts) / sum_wts

    # calculate how much extra shift is needed
    extra_z = spec.conv_blueshifts .- z_cbs_avg

    # loop over time
    p = Progress(Nt; enabled=show_progress)
    for t in 1:Nt
        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # trim all the bisector data
            @cusync @cuda threads=threads2 blocks=blocks2 trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                             depcontrast_gpu, lenall_gpu,
                                                                             bisall_gpu_loop, intall_gpu_loop,
                                                                             widall_gpu_loop, bisall_gpu,
                                                                             intall_gpu, widall_gpu, 
                                                                             bisall_mean, intall_mean, 
                                                                             widall_mean)

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

        # iterate progress meter
        next!(p)
    end

    # copy over flux
    @cusync flux_cpu .= Array(flux)

    # make sure nothing is still running on GPU
    # CUDA.synchronize()
    return nothing
end

function disk_sim_resolved_gpu(spec::SpecParams{T1}, disk::DiskParams{T1}, soldata::GPUSolarData{T2},
                               μ_bins::AA{T1,1}, gpu_allocs::GPUAllocsResolved{T2}, flux_cpu::AA{T1,3}; 
                               verbose::Bool=false, seed_rng::Bool=false,  
                               skip_times::BitVector=falses(disk.Nt),
                               show_progress::Bool=true) where {T1<:AF, T2<:AF}
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # parse out composite type
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux
    μ_bins_gpu = gpu_allocs.μ_bins_gpu

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

    # allocate destinations for interpolations
    @cusync begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
    end

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

    threads5 = (256, 4)
    blocks5 = (cld(CUDA.length(prof), prod(threads5)), cld(CUDA.length(μ_bins_gpu), prod(threads5)))

    # get weighted disk average cbs
    @cusync sum_wts = CUDA.sum(wts)
    @cusync z_cbs_avg = CUDA.sum(z_cbs .* wts) / sum_wts

    # calculate how much extra shift is needed
    extra_z = spec.conv_blueshifts .- z_cbs_avg

    # loop over time
    p = Progress(Nt; enabled=show_progress)
    for t in 1:Nt
        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
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
            @cusync @cuda threads=threads4 blocks=blocks4 line_profile_resolved_gpu!(prof, μ_bins_gpu, μs, wts, λs, allwavs, allints)

            # copy data from GPU to CPU
            @cusync @cuda threads=threads5 blocks=blocks5 apply_line_resolved!(t, prof, flux, sum_wts)
        end

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)

        # iterate progress meter
        next!(p)
    end

    # copy over flux
    # @cusync flux_cpu .= Array(flux)
    copyto!(flux_cpu, flux)

    # make sure nothing is still running on GPU
    # CUDA.synchronize()
    return nothing
end