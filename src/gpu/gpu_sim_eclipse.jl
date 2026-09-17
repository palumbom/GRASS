function disk_sim_eclipse_gpu(spec::SpecParams{T1}, disk::DiskParamsEclipse{T1}, 
                               soldata::GPUSolarData{T2}, gpu_allocs::GPUAllocsEclipse{T2},
                               flux_cpu::AA{T1,2}, obs_long::T1, obs_lat::T1, alt::T1, time_stamps::Vector{String}, wavelength,
                               ext_coeff, ext_toggle_gpu::Bool, spot_toggle_gpu::Bool, LD_type::String;
                               skip_times::BitVector=falses(disk.Nt),
                               interp_mu::Bool=false, pool_axes::Bool=false) where {T1<:AF, T2<:AF}

    # get dimensions for memory alloc
    Nt = disk.Nt
    Nλ = length(spec.lambdas)                           
                                
    # parse out GPU allocations
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux

    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints
    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs
    contrast = gpu_allocs.contrast

    # alias the input data from GPUSolarData
    lenall_gpu = soldata.len
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
    threads4 = (16,16)
    blocks4 = cld(CUDA.length(μs) * Nλ, prod(threads4))

    threads5 = 1024
    blocks5 = cld(CUDA.length(prof), prod(threads5))

    # allocate destinations for interpolations
    CUDA.@sync  begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
    end

    # allocate memory for means
    CUDA.@sync  begin
        bisall_mean = CUDA.zeros(CUDA.eltype(bisall_gpu_loop), 100, CUDA.size(bisall_gpu_loop, 3))
        intall_mean = CUDA.zeros(CUDA.eltype(intall_gpu_loop), 100, CUDA.size(intall_gpu_loop, 3))
        widall_mean = CUDA.zeros(CUDA.eltype(widall_gpu_loop), 100, CUDA.size(widall_gpu_loop, 3))
    end

    threads6 = (4, 16)
    blocks6 = cld(length(lenall_gpu) * 100, prod(threads6))

    CUDA.@sync  @cuda threads=threads6 blocks=blocks6 GRASS.time_average_bis!(lenall_gpu, bisall_mean, intall_mean,
                                                                    widall_mean, bisall_gpu, intall_gpu,
                                                                    widall_gpu)

    # interp_mu / pool_axes: time means of the TRIMMED arrays, recomputed per line inside the
    # epoch loop because trim_bisector_gpu! rewrites the _loop arrays per line, and their
    # axis-pooled copies. Allocated uninitialised: the kernels write every element, and a fill
    # here would be a kernel launch ahead of generate_tloop_gpu!, which changes the
    # granulation draw. With both flags off these alias the raw means and are never read.
    adjust_means = interp_mu | pool_axes
    if adjust_means
        bis_mean_trim = CuArray{CUDA.eltype(bisall_gpu)}(undef, 100, CUDA.size(bisall_gpu, 3))
        int_mean_trim = CuArray{CUDA.eltype(intall_gpu)}(undef, 100, CUDA.size(intall_gpu, 3))
        wid_mean_trim = CuArray{CUDA.eltype(widall_gpu)}(undef, 100, CUDA.size(widall_gpu, 3))
    else
        bis_mean_trim = bisall_mean
        int_mean_trim = intall_mean
        wid_mean_trim = widall_mean
    end
    if pool_axes
        bis_mean_src = CuArray{CUDA.eltype(bisall_gpu)}(undef, 100, CUDA.size(bisall_gpu, 3))
        int_mean_src = CuArray{CUDA.eltype(intall_gpu)}(undef, 100, CUDA.size(intall_gpu, 3))
        wid_mean_src = CuArray{CUDA.eltype(widall_gpu)}(undef, 100, CUDA.size(widall_gpu, 3))
    else
        bis_mean_src = bis_mean_trim
        int_mean_src = int_mean_trim
        wid_mean_src = wid_mean_trim
    end

    if ext_toggle_gpu == true
        ext_toggle_gpu = 1.0
    else 
        ext_toggle_gpu = 0.0
    end

    if spot_toggle_gpu == true
        spot_toggle_gpu = 1.0
    else 
        spot_toggle_gpu = 0.0
    end

    # loop over time
    for t in 1:Nt
        # sort out the system geometry
        calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, wavelength, LD_type, ext_toggle_gpu, ext_coeff, disk, gpu_allocs, spot_toggle_gpu)

        # get conv. blueshift and keys from input data
        get_keys_and_cbs_gpu!(gpu_allocs, soldata)

        if isone(t)
            # generate the random numbers on the gpu
            CUDA.@sync  @cuda threads=threads1 blocks=blocks1 GRASS.generate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
        end

        # interp_mu: bracketing tiles and limb-angle weights for this epoch's geometry, and
        # the interpolated z_cbs. MUST stay below generate_tloop_gpu!: that kernel draws
        # rand() on the device, the device RNG advances per kernel launch, and an extra launch
        # ahead of it changes the granulation realization.
        if interp_mu
            get_interp_keys_gpu!(gpu_allocs, soldata)
        end

        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            CUDA.@sync  @captured @cuda threads=threads1 blocks=blocks1 GRASS.iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            if ext_toggle_gpu == true
                # get weighted disk average cbs
                CUDA.@sync  sum_wts = CUDA.sum(dA .* ld[:,:,l] .* ext[:,:,l])
                CUDA.@sync  z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l] .* ext[:,:,l]) / sum_wts
            else
                # get weighted disk average cbs
                CUDA.@sync  sum_wts = CUDA.sum(dA .* ld[:,:,l])
                CUDA.@sync  z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l]) / sum_wts
            end

            # calculate how much extra shift is needed
            extra_z = spec.conv_blueshifts .- z_cbs_avg

            # trim all the bisector data
            CUDA.@sync  @cuda threads=threads2 blocks=blocks2 GRASS.trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                             depcontrast_gpu, lenall_gpu,
                                                                             bisall_gpu_loop, intall_gpu_loop,
                                                                             widall_gpu_loop, bisall_gpu,
                                                                             intall_gpu, widall_gpu)

            # interp_mu / pool_axes: time mean of the trimmed arrays, on their common depth
            # index, and the axis-pooled copy
            if adjust_means
                CUDA.@sync  @cuda threads=threads6 blocks=blocks6 GRASS.time_average_bis!(lenall_gpu, bis_mean_trim, int_mean_trim,
                                                                                wid_mean_trim, bisall_gpu_loop, intall_gpu_loop,
                                                                                widall_gpu_loop)
            end
            if pool_axes
                CUDA.@sync  @cuda threads=threads6 blocks=blocks6 GRASS.pool_axis_means_gpu!(bis_mean_src, int_mean_src, wid_mean_src,
                                                                                bis_mean_trim, int_mean_trim, wid_mean_trim,
                                                                                soldata.mu)
            end

            # assemble line shape on even int grid
            CUDA.@sync  @cuda threads=threads3 blocks=blocks3 fill_workspaces_2D_eclipse!(spec.lines[l], spec.variability[l], extra_z[l],
                                                                           tloop, dat_idx,
                                                                           z_rot, z_cbs, lenall_gpu,
                                                                           bisall_gpu_loop, intall_gpu_loop,
                                                                           widall_gpu_loop, allwavs, allints, contrast,
                                                                           interp_mu, pool_axes, gpu_allocs.dat_idx_lo, gpu_allocs.dat_idx_hi,
                                                                           gpu_allocs.dat_wt, bis_mean_src, int_mean_src, wid_mean_src,
                                                                           bis_mean_trim, int_mean_trim, wid_mean_trim)
            
            # do the line synthesis, interp back onto wavelength grid
            CUDA.@sync  @cuda threads=threads4 blocks=blocks4 GRASS.line_profile_gpu!(l, prof, μs, ld, dA, ext, λs, allwavs, allints, ext_toggle_gpu)

            # copy data from GPU to CPU
            CUDA.@sync  @cuda threads=threads5 blocks=blocks5 GRASS.apply_line!(t, prof, flux, sum_wts)
        end

        # iterate tloop
        CUDA.@sync  @captured @cuda threads=threads1 blocks=blocks1 GRASS.iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # copy over flux; skipped epochs are zero, as in the CPU eclipse path
    CUDA.@sync  flux_cpu .= Array(flux)
    flux_cpu[:, skip_times] .= zero(eltype(flux_cpu))

    # make sure nothing is still running on GPU
    CUDA.synchronize()
    return
end

function disk_sim_eclipse_gpu(spec::SpecParams{T1}, disk::DiskParamsEclipse{T1}, 
                               soldata::GPUSolarData{T2}, gpu_allocs::GPUAllocsEclipse{T2},
                               flux_cpu::AA{T1,2}, obs_long::T1, obs_lat::T1, alt::T1, time_stamps::Vector{String}, wavelength,
                               ext_coeff, CB1, CB2, CB3; skip_times::BitVector=falses(disk.Nt),
                               data_cbs::Bool=true,
                               static_bisector::Bool=false,
                               interp_mu::Bool=false, pool_axes::Bool=false) where {T1<:AF, T2<:AF}

    # get dimensions for memory alloc
    Nt = disk.Nt
    Nλ = length(spec.lambdas)                           
                                
    # parse out GPU allocations
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux

    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints
    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    contrast = gpu_allocs.contrast
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs

    # alias the input data from GPUSolarData
    lenall_gpu = soldata.len
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
    threads4 = (16,16)
    blocks4 = cld(CUDA.length(μs) * Nλ, prod(threads4))

    threads5 = 1024
    blocks5 = cld(CUDA.length(prof), prod(threads5))

    # allocate destinations for interpolations
    CUDA.@sync  begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
    end

    # allocate memory for means
    CUDA.@sync  begin
        bisall_mean = CUDA.zeros(CUDA.eltype(bisall_gpu_loop), 100, CUDA.size(bisall_gpu_loop, 3))
        intall_mean = CUDA.zeros(CUDA.eltype(intall_gpu_loop), 100, CUDA.size(intall_gpu_loop, 3))
        widall_mean = CUDA.zeros(CUDA.eltype(widall_gpu_loop), 100, CUDA.size(widall_gpu_loop, 3))
    end

    threads6 = (4, 16)
    blocks6 = cld(length(lenall_gpu) * 100, prod(threads6))

    CUDA.@sync  @cuda threads=threads6 blocks=blocks6 GRASS.time_average_bis!(lenall_gpu, bisall_mean, intall_mean,
                                                                    widall_mean, bisall_gpu, intall_gpu,
                                                                    widall_gpu)

    # interp_mu / pool_axes: time means of the TRIMMED arrays, recomputed per line inside the
    # epoch loop because trim_bisector_gpu! rewrites the _loop arrays per line, and their
    # axis-pooled copies. Allocated uninitialised: the kernels write every element, and a fill
    # here would be a kernel launch ahead of generate_tloop_gpu!, which changes the
    # granulation draw. With both flags off these alias the raw means and are never read.
    adjust_means = interp_mu | pool_axes
    if adjust_means
        bis_mean_trim = CuArray{CUDA.eltype(bisall_gpu)}(undef, 100, CUDA.size(bisall_gpu, 3))
        int_mean_trim = CuArray{CUDA.eltype(intall_gpu)}(undef, 100, CUDA.size(intall_gpu, 3))
        wid_mean_trim = CuArray{CUDA.eltype(widall_gpu)}(undef, 100, CUDA.size(widall_gpu, 3))
    else
        bis_mean_trim = bisall_mean
        int_mean_trim = intall_mean
        wid_mean_trim = widall_mean
    end
    if pool_axes
        bis_mean_src = CuArray{CUDA.eltype(bisall_gpu)}(undef, 100, CUDA.size(bisall_gpu, 3))
        int_mean_src = CuArray{CUDA.eltype(intall_gpu)}(undef, 100, CUDA.size(intall_gpu, 3))
        wid_mean_src = CuArray{CUDA.eltype(widall_gpu)}(undef, 100, CUDA.size(widall_gpu, 3))
    else
        bis_mean_src = bis_mean_trim
        int_mean_src = int_mean_trim
        wid_mean_src = wid_mean_trim
    end

    # static_bisector keeps the line asymmetry but removes its time evolution, by giving every
    # epoch of each tile the time-averaged profile that time_average_bis! just computed. That
    # separates the two things `variability = true` otherwise turns on together: the asymmetry
    # and the granulation jitter. Placed before the epoch loop, so it costs one launch and
    # cannot perturb the per-epoch RNG sequence.
    if static_bisector
        CUDA.@sync  @cuda threads=threads6 blocks=blocks6 GRASS.broadcast_mean_bis!(lenall_gpu,
                                                                    bisall_gpu, intall_gpu, widall_gpu,
                                                                    bisall_mean, intall_mean, widall_mean)
    end

    ext_toggle_gpu = 1.0
    # loop over time
    for t in 1:Nt
        calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, wavelength, ext_coeff, disk, gpu_allocs, CB1, CB2, CB3)

        # get conv. blueshift and keys from input data
        get_keys_and_cbs_gpu!(gpu_allocs, soldata)

        if isone(t)
            # generate the random numbers on the gpu
            CUDA.@sync  @cuda threads=threads1 blocks=blocks1 GRASS.generate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
        end

        # interp_mu: bracketing tiles and limb-angle weights for this epoch's geometry, and
        # the interpolated z_cbs (which the data_cbs zeroing below then overrides). MUST stay
        # below generate_tloop_gpu!, for the reason given at the data_cbs block.
        if interp_mu
            get_interp_keys_gpu!(gpu_allocs, soldata)
        end

        # data_cbs = false drops the data-driven convective blueshift while keeping the
        # bisector. Zero the array rather than gating the kernel term: the kernel applies one
        # flag to both z_cbs and extra_z, and extra_z carries spec.conv_blueshifts, so gating
        # there would silently discard a user-supplied blueshift. With z_cbs zeroed, z_cbs_avg
        # falls out to zero and extra_z reduces to exactly spec.conv_blueshifts.
        #
        # MUST stay below generate_tloop_gpu!. That kernel draws rand() on the device, and the
        # device RNG advances per kernel launch, so an extra launch ahead of it changes the
        # granulation realization. Placing this above it shifted Model III by 0.3 m/s while
        # leaving Models I and II untouched, which is how the ordering requirement was found.
        if !data_cbs
            CUDA.@sync  z_cbs .= zero(eltype(z_cbs))
        end

        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            CUDA.@sync  @captured @cuda threads=threads1 blocks=blocks1 GRASS.iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            if ext_toggle_gpu == true
                # get weighted disk average cbs
                CUDA.@sync  sum_wts = CUDA.sum(dA .* ld[:,:,l] .* ext[:,:,l])
                CUDA.@sync  z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l] .* ext[:,:,l]) / sum_wts
            else
                # get weighted disk average cbs
                CUDA.@sync  sum_wts = CUDA.sum(dA .* ld[:,:,l])
                CUDA.@sync  z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l]) / sum_wts
            end

            # calculate how much extra shift is needed
            extra_z = spec.conv_blueshifts .- z_cbs_avg

            # trim all the bisector data
            CUDA.@sync  @cuda threads=threads2 blocks=blocks2 GRASS.trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                             depcontrast_gpu, lenall_gpu,
                                                                             bisall_gpu_loop, intall_gpu_loop,
                                                                             widall_gpu_loop, bisall_gpu,
                                                                             intall_gpu, widall_gpu)

            # interp_mu / pool_axes: time mean of the trimmed arrays, on their common depth
            # index, and the axis-pooled copy
            if adjust_means
                CUDA.@sync  @cuda threads=threads6 blocks=blocks6 GRASS.time_average_bis!(lenall_gpu, bis_mean_trim, int_mean_trim,
                                                                                wid_mean_trim, bisall_gpu_loop, intall_gpu_loop,
                                                                                widall_gpu_loop)
            end
            if pool_axes
                CUDA.@sync  @cuda threads=threads6 blocks=blocks6 GRASS.pool_axis_means_gpu!(bis_mean_src, int_mean_src, wid_mean_src,
                                                                                bis_mean_trim, int_mean_trim, wid_mean_trim,
                                                                                soldata.mu)
            end

            # assemble line shape on even int grid. `variability` gates the bisector (line
            # asymmetry and its time evolution) in trim_bisector_gpu! above, and here it gates
            # z_cbs and the extra_z that re-references it. The data_cbs switch acts on the
            # z_cbs array above, not on this flag, so spec.conv_blueshifts survives it.
            CUDA.@sync  @cuda threads=threads3 blocks=blocks3 fill_workspaces_2D_eclipse!(spec.lines[l], spec.variability[l], extra_z[l],
                                                                           tloop, dat_idx,
                                                                           z_rot, z_cbs, lenall_gpu,
                                                                           bisall_gpu_loop, intall_gpu_loop,
                                                                           widall_gpu_loop, allwavs, allints, contrast,
                                                                           interp_mu, pool_axes, gpu_allocs.dat_idx_lo, gpu_allocs.dat_idx_hi,
                                                                           gpu_allocs.dat_wt, bis_mean_src, int_mean_src, wid_mean_src,
                                                                           bis_mean_trim, int_mean_trim, wid_mean_trim)
            
            # do the line synthesis, interp back onto wavelength grid
            CUDA.@sync  @cuda threads=threads4 blocks=blocks4 GRASS.line_profile_gpu!(l, prof, μs, ld, dA, ext, λs, allwavs, allints, ext_toggle_gpu)

            # copy data from GPU to CPU
            CUDA.@sync  @cuda threads=threads5 blocks=blocks5 GRASS.apply_line!(t, prof, flux, sum_wts)
        end

        # iterate tloop
        CUDA.@sync  @captured @cuda threads=threads1 blocks=blocks1 GRASS.iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # copy over flux; skipped epochs are zero, as in the CPU eclipse path
    CUDA.@sync  flux_cpu .= Array(flux)
    flux_cpu[:, skip_times] .= zero(eltype(flux_cpu))

    # make sure nothing is still running on GPU
    CUDA.synchronize()
    return
end
