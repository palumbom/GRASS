function iterate_tloop_gpu(tloop, data_inds, lenall, grid)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            # find position on disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)

            # move to next iter if off disk
            if r2 > 1.0
                continue
            end

            # iterate tloop
            @inbounds tloop[i,j] = tloop[i,j] + 1

            # check that tloop didn't overshoot the data
            ntimes = lenall[data_inds[i,j]]
            if tloop[i,j] >= ntimes
                @inbounds tloop[i,j] = 1
            end
        end
    end
    return nothing
end

function initialize_arrays_for_gpu(data_inds, tloop, norm_term, grid, disc_mu, disc_ax, lenall, u1, u2)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            # find position on disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)

            # move to next iter if off disk
            if r2 > 1.0
                continue
            end

            # find the nearest mu ind and ax code
            mu = calc_mu(r2)
            nn_mu_ind = searchsortednearest_gpu(disc_mu, mu)
            nn_ax_code = find_nearest_ax_gpu(x, y)

            # find the correct data index
            @inbounds data_inds[i,j] = find_data_index_gpu(nn_mu_ind, nn_ax_code)

            # iterate tloop
            @inbounds tloop[i,j] = tloop[i,j] + 1

            # check that tloop didn't overshoot the data
            ntimes = lenall[data_inds[i,j]]
            if tloop[i,j] >= ntimes
                @inbounds tloop[i,j] = 1
            end

            # calculate the normalization
            @inbounds norm_term[i,j] = calc_norm_term(mu, CUDA.length(grid), u1, u2)
        end
    end
    return nothing
end

function disk_sim_gpu(spec::SpecParams, disk::DiskParams, outspec::AA{T,2}) where T<:Float64
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # make array to copy GPU results into
    starmap_cpu = zeros(N, N, Nλ)

    # sort the input data for use on GPU
    sorted_data = sort_data_for_gpu(spec.soldata)
    disc_mu = sorted_data[1]
    disc_ax = sorted_data[2]
    lenall_cpu = sorted_data[3]
    wavall_cpu = sorted_data[4]
    bisall_cpu = sorted_data[5]
    widall_cpu = sorted_data[6]
    depall_cpu = sorted_data[7]

    # move input data to gpu
    disc_mu = CuArray(disc_mu)
    disc_ax = CuArray(disc_ax)
    lenall_gpu = CuArray(lenall_cpu)
    wavall_gpu = CuArray(wavall_cpu)
    bisall_gpu = CuArray(bisall_cpu)
    widall_gpu = CuArray(widall_cpu)
    depall_gpu = CuArray(depall_cpu)

    # allocate arrays for fresh copy of input data to copy to each loop
    wavall_gpu_loop = CUDA.copy(wavall_gpu)
    bisall_gpu_loop = CUDA.copy(bisall_gpu)
    widall_gpu_loop = CUDA.copy(widall_gpu)
    depall_gpu_loop = CUDA.copy(depall_gpu)

    # allocate memory for synthesis on the GPU
    starmap = CUDA.ones(Float64, N, N, Nλ)
    lwavgrid = CUDA.zeros(Float64, N, N, 100)
    rwavgrid = CUDA.zeros(Float64, N, N, 100)
    allwavs = CUDA.zeros(Float64, N, N, 200)
    allints = CUDA.zeros(Float64, N, N, 200)

    # get random starting indices and move it to GPU
    tloop = rand(1:maximum(lenall_cpu), (N, N))
    tloop = CuArray(tloop)

    # move other data to the gpu
    grid = CuArray(GRASS.make_grid(N))
    lambdas = CuArray(spec.lambdas)

    # allocate memory for input data indices and normalization terms
    data_inds = CuArray{Int64,2}(undef, N, N)
    norm_term = CuArray{Float64,2}(undef, N, N)

    # set number of threads and blocks for N*N matrix gpu functions
    threads2 = (16, 16)
    blocks2 = cld(N^2, prod(threads2))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads3 = (6,6,6)
    blocks3 = cld(N^2 * Nλ, prod(threads3))

    # initialize values for data_inds, tloop, and norm_term
    CUDA.@sync @cuda threads=threads2 blocks=blocks2 initialize_arrays_for_gpu(data_inds, tloop, norm_term, grid, disc_mu,
                                                                               disc_ax, lenall_gpu, disk.u1, disk.u2)

    # loop over time
    for t in 1:Nt
        # initialize starmap with fresh copy of weights
        starmap .= CUDA.copy(norm_term)

        # copy the clean, untrimmed input data to workspace each time iteration
        # TODO does this need to be moved down a loop?
        CUDA.@sync begin
            wavall_gpu_loop .= CUDA.copy(wavall_gpu)
            bisall_gpu_loop .= CUDA.copy(bisall_gpu)
            widall_gpu_loop .= CUDA.copy(widall_gpu)
            depall_gpu_loop .= CUDA.copy(depall_gpu)
        end

        # loop over lines to synthesize
        for l in 1:length(spec.lines)
            # pre-trim the data
            for n in 1:length(lenall_cpu)
                CUDA.@sync begin
                    # send slices to gpu
                    wavall_gpu_loop_slice = CUDA.view(wavall_gpu_loop, :, 1:lenall_cpu[n], n)
                    bisall_gpu_loop_slice = CUDA.view(bisall_gpu_loop, :, 1:lenall_cpu[n], n)
                    widall_gpu_loop_slice = CUDA.view(widall_gpu_loop, :, 1:lenall_cpu[n], n)
                    depall_gpu_loop_slice = CUDA.view(depall_gpu_loop, :, 1:lenall_cpu[n], n)

                    # do the trim
                    threads1 = 100
                    blocks1 = cld(lenall_cpu[n] * 100, threads1)
                    @cuda threads=threads1 blocks=blocks1 trim_bisector_chop_gpu!(spec.depths[l],
                                                                                  wavall_gpu_loop_slice,
                                                                                  bisall_gpu_loop_slice,
                                                                                  depall_gpu_loop_slice,
                                                                                  widall_gpu_loop_slice,
                                                                                  NaN)
                end
            end

            # do the line synthesis
            CUDA.@sync @cuda threads=threads3 blocks=blocks3 line_profile_gpu!(starmap, tloop, spec.lines[l],
                                                                               spec.depths[l], spec.conv_blueshifts[l],
                                                                               grid, lambdas, data_inds, lenall_gpu,
                                                                               wavall_gpu_loop, bisall_gpu_loop,
                                                                               widall_gpu_loop, depall_gpu_loop,
                                                                               lwavgrid, rwavgrid,allwavs, allints)

            # do array reduction and move data from GPU to CPU
            CUDA.@sync starmap_cpu .= Array(CUDA.sum(starmap, dims=(1,2)))
            outspec[:,t] .= view(starmap_cpu, 1, 1, :)

            # iterate tloop
            if t < Nt
                CUDA.@sync @cuda threads=threads2 blocks=blocks2 iterate_tloop_gpu(tloop, data_inds, lenall_gpu, grid)
            end
        end
    end
    return spec.lambdas, outspec
end
