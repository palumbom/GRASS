function disk_sim_gpu(spec::SpecParams{T1}, disk::DiskParams{T1}, soldata::GPUSolarData{T2},
                      gpu_allocs::GPUAllocs{T2}, outspec::AA{T1,2}; verbose::Bool=false,
                      seed_rng::Bool=false,  skip_times::BitVector=falses(disk.Nt)) where {T1<:AF, T2<:AF}
    # infer precision
    precision = T2

    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # parse out composite type
    λs = gpu_allocs.λs
    μs = gpu_allocs.μs
    wts = gpu_allocs.wts
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs

    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    starmap = gpu_allocs.starmap
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

    # set number of threads and blocks for N*N matrix gpu functions
    threads1 = (16,16)
    blocks1 = cld(N^2, prod(threads1))

    # set number of threads and blocks for trimming functions
    threads2 = (4,4,16)
    blocks2 = cld(length(lenall_gpu) * maximum(lenall_gpu) * 100, prod(threads2))

    # set number of threads and blocks for N*N*100 matrix gpu functions
    threads3 = (4,4,16)
    blocks3 = cld(N^2 * 100, prod(threads3))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads4 = (3,3,42)
    blocks4 = cld(N^2 * Nλ, prod(threads4))

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        bisall_gpu_loop = CUDA.zeros(precision, CUDA.size(bisall_gpu))
        intall_gpu_loop = CUDA.zeros(precision, CUDA.size(intall_gpu))
        widall_gpu_loop = CUDA.zeros(precision, CUDA.size(widall_gpu))
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
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, μs, dat_idx, lenall_gpu)
            continue
        end

        # initialize starmap with fresh copy of weights
        @cusync starmap .= wts

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            @cusync begin
                CUDA.copyto!(bisall_gpu_loop, bisall_gpu)
                CUDA.copyto!(intall_gpu_loop, intall_gpu)
                CUDA.copyto!(widall_gpu_loop, widall_gpu)
            end

            # trim all the bisector data
            @cusync @captured @cuda threads=threads2 blocks=blocks2 trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                                       depcontrast_gpu, lenall_gpu,
                                                                                       bisall_gpu_loop, intall_gpu_loop,
                                                                                       widall_gpu_loop, bisall_gpu,
                                                                                       intall_gpu, widall_gpu)

            # assemble line shape on even int grid
            @cusync @captured @cuda threads=threads3 blocks=blocks3 fill_workspaces!(spec.lines[l], spec.variability[l],
                                                                                     extra_z[l], μs, tloop, dat_idx,
                                                                                     z_rot, z_cbs, bisall_gpu_loop,
                                                                                     intall_gpu_loop, widall_gpu_loop,
                                                                                     allwavs, allints)

            # do the line synthesis, interp back onto wavelength grid
            @cusync @captured @cuda threads=threads4 blocks=blocks4 line_profile_gpu!(starmap, μs, λs, allwavs, allints)
        end

        # do array reduction and move data from GPU to CPU
        @cusync @inbounds outspec[:,t] .*= dropdims(Array(CUDA.sum(starmap, dims=(1,2))), dims=(1,2))

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, μs, dat_idx, lenall_gpu)
    end

    # ensure normalization
    outspec ./= sum_wts
    CUDA.synchronize()
    return nothing
end
