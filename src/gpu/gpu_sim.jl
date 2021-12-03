function iterate_tloop_gpu(tloop, data_inds, lenall, grid)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            # find position on disk and move to next iter if off disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
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

function initialize_arrays_for_gpu(data_inds, norm_terms, rot_shifts,
                                   grid, disc_mu, disc_ax, u1, u2,
                                   polex, poley, polez)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            # find position on disk and move to next iter if off disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
            if r2 > 1.0
                continue
            end

            # find the nearest mu ind and ax code
            mu = calc_mu(r2)
            nn_mu_ind = searchsortednearest_gpu(disc_mu, mu)
            nn_ax_code = find_nearest_ax_gpu(x, y)

            # find the correct data index
            @inbounds data_inds[i,j] = find_data_index_gpu(nn_mu_ind, nn_ax_code)

            # calculate the normalization
            @inbounds norm_terms[i,j] = calc_norm_term(mu, CUDA.length(grid), u1, u2)

            # calculate the rotational doppler shift
            @inbounds rot_shifts[i,j] = patch_velocity_los_gpu(x, y, 1.0, polex, poley, polez)
        end
    end
    return nothing
end

function disk_sim_gpu(spec::SpecParams, disk::DiskParams, outspec::AA{T,2}; skip_times::BitVector=BitVector(zeros(disk.Nt))) where T<:Float64
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # get pole component vectors
    polex, poley, polez = disk.pole

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
    CUDA.@sync begin
        disc_mu = CuArray(disc_mu)
        disc_ax = CuArray(disc_ax)
        lenall_gpu = CuArray(lenall_cpu)
        wavall_gpu = CuArray(wavall_cpu)
        bisall_gpu = CuArray(bisall_cpu)
        widall_gpu = CuArray(widall_cpu)
        depall_gpu = CuArray(depall_cpu)
    end

    # allocate arrays for fresh copy of input data to copy to each loop
    CUDA.@sync begin
        wavall_gpu_loop = CUDA.copy(wavall_gpu)
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
        depall_gpu_loop = CUDA.copy(depall_gpu)
    end

    # allocate memory for synthesis on the GPU
    starmap_cpu = zeros(N, N, Nλ)
    CUDA.@sync begin
        starmap = CUDA.ones(Float64, N, N, Nλ)
        lwavgrid = CUDA.zeros(Float64, N, N, 100)
        rwavgrid = CUDA.zeros(Float64, N, N, 100)
        allwavs = CUDA.zeros(Float64, N, N, 200)
        allints = CUDA.zeros(Float64, N, N, 200)
    end

    # get random starting indices and move it to GPU
    tloop = rand(1:maximum(lenall_cpu), (N, N))
    # tloop = ones(Int64,N,N)
    CUDA.@sync tloop = CuArray(tloop)

    # move other data to the gpu
    CUDA.@sync begin
        grid = CuArray(GRASS.make_grid(N))
        lambdas = CuArray(spec.lambdas)
    end

    # allocate memory for input data indices and normalization terms
    CUDA.@sync begin
        data_inds = CuArray{Int64,2}(undef, N, N)
        norm_terms = CuArray{Float64,2}(undef, N, N)
        rot_shifts = CuArray{Float64,2}(undef, N, N)
        λΔDs = CuArray{Float64,2}(undef, N, N)
    end

    # set number of threads and blocks for N*N matrix gpu functions
    threads2 = (16, 16)
    blocks2 = cld(N^2, prod(threads2))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads3 = (6,6,6)
    blocks3 = cld(N^2 * Nλ, prod(threads3))

    # initialize values for data_inds, tloop, and norm_terms
    CUDA.@sync @cuda threads=threads2 blocks=blocks2 initialize_arrays_for_gpu(data_inds, norm_terms, rot_shifts,
                                                                               grid, disc_mu, disc_ax, disk.u1,
                                                                               disk.u2, polex, poley, polez)

    # check that random indices for time index don't exceed dataset length
    CUDA.@sync @cuda threads=threads2 blocks=blocks2 iterate_tloop_gpu(tloop, data_inds, lenall_gpu, grid)

    # loop over time
    for t in 1:Nt
        # don't do all this work if skip_times is true
        if skip_times[t]
            continue
        end

        # initialize starmap with fresh copy of weights
        # CUDA.@sync starmap .= CUDA.copy(norm_terms)
        CUDA.@sync starmap .= 1.0

        # copy the clean, untrimmed input data to workspace each time iteration
        CUDA.@sync begin
            wavall_gpu_loop .= CUDA.copy(wavall_gpu)
            bisall_gpu_loop .= CUDA.copy(bisall_gpu)
            widall_gpu_loop .= CUDA.copy(widall_gpu)
            depall_gpu_loop .= CUDA.copy(depall_gpu)
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # pre-trim the data, loop over all disk positions
            for n in eachindex(lenall_cpu)
                CUDA.@sync begin
                    # get correct index for position of input data
                    wavall_gpu_out = CUDA.view(wavall_gpu_loop, :, 1:lenall_cpu[n], n) #.* spec.variability[l]
                    bisall_gpu_out = CUDA.view(bisall_gpu_loop, :, 1:lenall_cpu[n], n)
                    widall_gpu_out = CUDA.view(widall_gpu_loop, :, 1:lenall_cpu[n], n)
                    depall_gpu_out = CUDA.view(depall_gpu_loop, :, 1:lenall_cpu[n], n)

                    wavall_gpu_in = CUDA.view(wavall_gpu, :, 1:lenall_cpu[n], n) .* spec.variability[l]
                    bisall_gpu_in = CUDA.view(bisall_gpu, :, 1:lenall_cpu[n], n)
                    widall_gpu_in = CUDA.view(widall_gpu, :, 1:lenall_cpu[n], n)
                    depall_gpu_in = CUDA.view(depall_gpu, :, 1:lenall_cpu[n], n)
                end

                # do the trim
                threads1 = (16,16)
                blocks1 = cld(lenall_cpu[n] * 100, prod(threads1))
                CUDA.@sync @cuda threads=threads1 blocks=blocks1 trim_bisector_chop_gpu(spec.depths[l],
                                                                                        wavall_gpu_out, bisall_gpu_out,
                                                                                        depall_gpu_out, widall_gpu_out,
                                                                                        wavall_gpu_in, bisall_gpu_in,
                                                                                        depall_gpu_in, widall_gpu_in,
                                                                                        NaN)
            end

            # fill workspace arrays
            # TODO: compare kernel launch cost to compute cost
            # TODO: how many registers are needed?
            # change here
            threads4 = (6,6,6)
            blocks4 = cld(N^2 * 100, prod(threads4))
            CUDA.@sync @cuda threads=threads4 blocks=blocks4 fill_workspace_arrays!(spec.lines[l], spec.depths[l],
                                                                                    spec.conv_blueshifts[l], grid, tloop,
                                                                                    data_inds, rot_shifts, λΔDs,
                                                                                    lenall_gpu, wavall_gpu_loop,
                                                                                    bisall_gpu_loop, widall_gpu_loop,
                                                                                    depall_gpu_loop, lwavgrid,
                                                                                    rwavgrid, allwavs, allints)
            CUDA.@sync @cuda threads=threads4 blocks=blocks4 concatenate_workspace_arrays!(spec.lines[l], spec.depths[l],
                                                                                           spec.conv_blueshifts[l], grid, tloop,
                                                                                           data_inds, rot_shifts, λΔDs,
                                                                                           lenall_gpu, wavall_gpu_loop,
                                                                                           bisall_gpu_loop, widall_gpu_loop,
                                                                                           depall_gpu_loop, lwavgrid,
                                                                                           rwavgrid, allwavs, allints)

            # do the line synthesis
            CUDA.@sync @cuda threads=threads3 blocks=blocks3 line_profile_gpu!(starmap, tloop, spec.lines[l],
                                                                               spec.depths[l], spec.conv_blueshifts[l],
                                                                               grid, lambdas, data_inds, rot_shifts, λΔDs,
                                                                               allwavs, allints)

            # do array reduction and move data from GPU to CPU
            CUDA.@sync starmap_cpu .= Array(starmap .* norm_terms)
            outspec[:,t] .= dropdims(sum(starmap_cpu, dims=(1,2)), dims=(1,2))
            # CUDA.@sync outspec[:,t] .= view(Array(CUDA.sum(starmap, dims=(1,2))), 1, 1, :)

            # iterate tloop
            if t < Nt
                CUDA.@sync @cuda threads=threads2 blocks=blocks2 iterate_tloop_gpu(tloop, data_inds, lenall_gpu, grid)
            end
        end
    end
    println("derp2")
    return spec.lambdas, outspec
end
