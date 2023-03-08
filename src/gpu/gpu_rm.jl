function disk_sim_rm_gpu(spec::SpecParams{T}, disk::DiskParams{T}, planet::Planet,
                         soldata::SolarData{T}, gpu_allocs::GPUAllocs, outspec::AA{T,2};
                         verbose::Bool=false, precision::DataType=Float64,
                         skip_times::BitVector=BitVector(zeros(disk.Nt))) where T<:AF
    # set single or double precision
    prec = precision

    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # get pole component vectors and limb darkening parameters
    polex, poley, polez = disk.pole
    u1 = disk.u1
    u2 = disk.u2

    # parse out composite type
    grid = gpu_allocs.grid
    lambdas = gpu_allocs.lambdas
    tloop = gpu_allocs.tloop
    data_inds = gpu_allocs.data_inds
    norm_terms = gpu_allocs.norm_terms
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs
    starmap = gpu_allocs.starmap
    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints

    # sort the input data for use on GPU
    sorted_data = sort_data_for_gpu(soldata)
    disc_mu_cpu = sorted_data[1]
    disc_ax_cpu = sorted_data[2]
    lenall_cpu = sorted_data[3]
    cbsall_cpu = sorted_data[4]
    bisall_cpu = sorted_data[5]
    intall_cpu = sorted_data[6]
    widall_cpu = sorted_data[7]

    # set number of threads and blocks for N*N matrix gpu functions
    threads1 = (16,16)
    blocks1 = cld(N^2, prod(threads1))

    # set number of threads and blocks for trimming functions
    threads2 = (4,4,16)
    blocks2 = cld(length(lenall_cpu) * maximum(lenall_cpu) * 100, prod(threads2))

    # set number of threads and blocks for N*N*100 matrix gpu functions
    threads3 = (4,4,16)
    blocks3 = cld(N^2 * 100, prod(threads3))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads4 = (3,3,42)
    blocks4 = cld(N^2 * Nλ, prod(threads4))

    # move input data to gpu
    @cusync begin
        disc_mu_gpu = CuArray{prec}(disc_mu_cpu)
        disc_ax_gpu = CuArray{Int32}(disc_ax_cpu)
        lenall_gpu = CuArray{Int32}(lenall_cpu)
        cbsall_gpu = CuArray{prec}(cbsall_cpu)
        bisall_gpu = CuArray{prec}(bisall_cpu)
        intall_gpu = CuArray{prec}(intall_cpu)
        widall_gpu = CuArray{prec}(widall_cpu)
    end

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
    end

    # # get launch parameters
    # kernel = @cuda launch=false initialize_arrays_for_gpu(tloop, data_inds, norm_terms, z_rot, z_cbs,
    #                                                       grid, disc_mu_gpu, disc_ax_gpu, lenall_gpu,
    #                                                       cbsall_gpu, u1, u2, polex, poley, polez)

    # config = launch_configuration(kernel.fun)
    # threads1 = min(N^2, config.threads)
    # blocks1 = cld(N^2, threads1)

    # @show threads1
    # @show blocks1

    # kernel(tloop, data_inds, norm_terms, z_rot, z_cbs, grid,
    #        disc_mu_gpu, disc_ax_gpu, lenall_gpu, cbsall_gpu,
    #        u1, u2, polex, poley, polez; threads=threads1, blocks=blocks1)


    # initialize values for data_inds, tloop, dop_shifts, and norm_terms
    @cusync @cuda threads=threads1 blocks=blocks1 initialize_arrays_for_gpu(tloop, data_inds, norm_terms,
                                                                            z_rot, z_cbs, grid, disc_mu_gpu,
                                                                            disc_ax_gpu, lenall_gpu, cbsall_gpu,
                                                                            u1, u2, polex, poley, polez)

    # get weighted disk average cbs
    @cusync sum_norm_terms = CUDA.sum(norm_terms)
    @cusync z_cbs_avg = CUDA.sum(z_cbs .* norm_terms) / sum_norm_terms

    # calculate how much extra shift is needed
    extra_z = spec.conv_blueshifts .- z_cbs_avg

    # loop over time
    for t in 1:Nt
        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, data_inds, lenall_gpu, grid)
            continue
        end

        # initialize starmap with fresh copy of weights
        @cusync starmap .= norm_terms

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # trim all the bisector data
            @cusync @captured @cuda threads=threads2 blocks=blocks2 trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                                       lenall_gpu, bisall_gpu_loop,
                                                                                       intall_gpu_loop, bisall_gpu,
                                                                                       intall_gpu)

            # assemble line shape on even int grid
            @cusync @captured @cuda threads=threads3 blocks=blocks3 fill_workspaces!(spec.lines[l], extra_z[l], grid,
                                                                                     tloop, data_inds, z_rot, z_cbs,
                                                                                     bisall_gpu_loop, intall_gpu_loop,
                                                                                     widall_gpu, allwavs, allints)

            # do the line synthesis, interp back onto wavelength grid
            @cusync @captured @cuda threads=threads4 blocks=blocks4 line_profile_gpu!(starmap, grid, lambdas, allwavs, allints)
        end

        # do array reduction and move data from GPU to CPU
        @cusync @inbounds outspec[:,t] .*= dropdims(Array(CUDA.sum(starmap, dims=(1,2))), dims=(1,2))

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, data_inds, lenall_gpu, grid)
    end

    # ensure normalization
    outspec ./= sum_norm_terms
    CUDA.synchronize()
    return nothing
end

# TODO: allocate the most I'll need
# TODO: add up whole sun where i,j are normal pixels
# TODO: then calculate profiles for occulted cells and subtract off
function calc_norm_term_rm_gpu(star_map, norm_terms, rot_shifts, grid, xpos, ypos, r_planet, u1, u2, polex, poley, polez)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            @inbounds star_map[i,j,:] .= norm_terms[i,j]
        end
    end
    return nothing
end
