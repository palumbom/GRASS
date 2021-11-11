function line_profile_gpu!(mid, lambdas, prof, wavm, depm, widm, lwavgrid, rwavgrid, allwavs, allints)
    # set wavgrids to line center to start
    for i in 1:CUDA.length(lwavgrid)
        lwavgrid[i] = (mid - (0.5 * widm[i] - wavm[i]))
        rwavgrid[i] = (mid + (0.5 * widm[i] + wavm[i]))
    end
    rwavgrid[1] = lwavgrid[1] + 1e-3            # TODO: fix to deal with nodes

    # concatenate into one big array
    len = CUDA.length(rwavgrid)
    for i in 1:CUDA.length(rwavgrid)
        allwavs[i+len] = rwavgrid[i]
        allints[i+len] = depm[i]
        allwavs[i] = lwavgrid[CUDA.length(rwavgrid) - (i - 1)]
        allints[i] = depm[CUDA.length(rwavgrid) - (i - 1)]
    end

    # interpolate onto original lambda grid, extrapolate to continuum
    linear_interp_mult_gpu(prof, lambdas, allwavs, allints, 1.0)
    return nothing
end

function calc_data_ind_gpu(data_inds, grid, disc_mu, disc_ax)
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
        end
    end
    return nothing
end

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

function disk_sim(star_map, tloop, lines, depths, z_convs, grid,
                  lambdas, data_inds, lenall, wavall, bisall, widall,
                  depall, lwavgrid, rwavgrid, allwavs, allints)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y
    idz = threadIdx().z + blockDim().z * (blockIdx().z-1)
    sdz = blockDim().z * gridDim().z

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            for k in idz:sdz:CUDA.length(lambdas)
                # set intensity to zero and go to next iter if off grid
                x = grid[i]
                y = grid[j]
                r2 = calc_r2(x, y)
                if r2 > 1.0
                    @inbounds star_map[i,j,k] = 0.0
                    continue
                end

                # calculate the shifted center of the line
                z_rot = patch_velocity_los(x, y)
                λΔD = lines * (1.0 + z_rot) * (1.0 + z_convs)

                # skip all this work if far from line core
                if (lambdas[k] < (λΔD - 0.5)) | (lambdas[k] > (λΔD + 0.5))
                    continue
                end

                # slice out the correct views of the input data for position
                data_ind = data_inds[i,j]
                wavt = CUDA.view(wavall, :, tloop[i,j], data_ind)
                bist = CUDA.view(bisall, :, tloop[i,j], data_ind)
                widt = CUDA.view(widall, :, tloop[i,j], data_ind)
                dept = CUDA.view(depall, :, tloop[i,j], data_ind)

                # set wavgrids based on bisector + wid data
                for n in 1:CUDA.size(lwavgrid,3)
                    @inbounds lwavgrid[i,j,n] = (λΔD - (0.5 * widt[n] - wavt[n]))
                    @inbounds rwavgrid[i,j,n] = (λΔD + (0.5 * widt[n] + wavt[n]))
                end
                @inbounds rwavgrid[i,j,1] = lwavgrid[i,j,1] + 1e-3

                # concatenate wavgrids into one big array
                len = CUDA.size(rwavgrid,3)
                for n in 1:CUDA.size(lwavgrid,3)
                    @inbounds allwavs[i,j,n+len] = rwavgrid[i,j,n]
                    @inbounds allints[i,j,n+len] = dept[n]
                    @inbounds allwavs[i,j,n] = lwavgrid[i,j, CUDA.size(rwavgrid,3) - (n - 1)]
                    @inbounds allints[i,j,n] = dept[CUDA.size(rwavgrid,3) - (n - 1)]
                end

                # take view of arrays to pass to interpolater
                allwavs_ij = CUDA.view(allwavs, i, j, :)
                allints_ij = CUDA.view(allints, i, j, :)

                # do the interpolation
                if ((lambdas[k] < CUDA.first(allwavs_ij)) | (lambdas[k] > CUDA.last(allwavs_ij)))
                    factor = 1.0
                else
                    m = CUDA.searchsortedfirst(allwavs_ij, lambdas[k]) - 1
                    m0 = CUDA.clamp(m, CUDA.firstindex(allints_ij), CUDA.lastindex(allints_ij))
                    m1 = CUDA.clamp(m+1, CUDA.firstindex(allints_ij), CUDA.lastindex(allints_ij))
                    factor = (allints_ij[m0] + (allints_ij[m1] - allints_ij[m0]) * (lambdas[k] - allwavs_ij[m0]) / (allwavs_ij[m1] - allwavs_ij[m0]))
                end
                @inbounds star_map[i,j,k] = star_map[i,j,k] * factor
            end
        end
    end
    return nothing
end
