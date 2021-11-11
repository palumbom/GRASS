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
            tloop[i,j] = tloop[i,j] + 1

            # check that tloop didn't overshoot the data
            ntimes = lenall[data_inds[i,j]]
            if tloop[i,j] >= ntimes
                @inbounds tloop[i,j] = 1
            end
        end
    end
    return nothing
end

function calc_norm_term_gpu(star_map, grid, u1, u2)
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
                @inbounds star_map[i,j,:] .= 0.0
                continue
            end

            # calculate the normalization
            mu = calc_mu(r2)
            @inbounds star_map[i,j,:] .= calc_norm_term(mu, CUDA.length(grid), u1, u2)
        end
    end
    return nothing
end

function disk_sim_gpu(spec, disk, outspec)
    # get random starting indices
    tloop = rand(1:maximum(lenall_cpu), (N, N))

    # allocate memory for input data indices
    data_inds = CuArray{Int64,2}(undef, N, N)

    # before starting loops, get map of data indices
    threads = (16, 16)
    blocks = cld(N^2, prod(threads))
    CUDA.@sync @cuda threads=threads blocks=blocks calc_data_ind_gpu(data_inds, grid, disc_mu, disc_ax)

    # make tloop starting index doesnt exceed number of epochs for position
    threads = (16, 16)
    blocks = cld(N^2, prod(threads))
    CUDA.@sync @cuda threads=threads blocks=blocks iterate_tloop_gpu(tloop, data_inds, lenall_gpu, grid)

    # loop over time
    for t in 1:Nt
        # initialize star_map with normalization term
        threads1 = (16, 16)
        blocks1 = cld(N^2, prod(threads1))
        CUDA.@sync @cuda threads=threads1 blocks=blocks1 calc_norm_term_gpu(starmap, grid, 0.4, 0.26)

        # copy the clean, untrimmed input data to workspace each time iteration
        # TODO does this need to be moved down a loop?
        CUDA.@sync begin
            wavall_gpu_loop .= CUDA.copy(wavall_gpu)
            bisall_gpu_loop .= CUDA.copy(bisall_gpu)
            widall_gpu_loop .= CUDA.copy(widall_gpu)
            depall_gpu_loop .= CUDA.copy(depall_gpu)
        end

        # loop over lines to synthesize
        for l in 1:length(lines)
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
            threads1 = (6,6,6)
            blocks1 = cld(N^2 * Nλ, prod(threads1))
            CUDA.@sync @cuda threads=threads1 blocks=blocks1 line_profile_gpu!(starmap, tloop, spec.lines[l],
                                                                               spec.depths[l], spec.conv_blueshifts[l],
                                                                               grid, lambdas, data_inds, lenall_gpu,
                                                                               wavall_gpu_loop, bisall_gpu_loop,
                                                                               widall_gpu_loop, depall_gpu_loop,
                                                                               lwavgrid, rwavgrid,allwavs, allints)

            # do array reduction and move data from GPU to CPU
            CUDA.@sync starmap_cpu_copy = Array(CUDA.sum(starmap, dims=(1,2)))
            outspec[:,t] .= view(starmap_cpu_copy, 1, 1, :)

            # iterate tloop
            if t < Nt
                threads1 = (16, 16)
                blocks1 = cld(N^2, prod(threads1))
                CUDA.@sync @cuda threads=threads1 blocks=blocks1 GRASS.iterate_tloop_gpu(tloop, data_inds, lenall_gpu, grid)
            end
        end
    end
    return spec.lambdas, outspec
end
