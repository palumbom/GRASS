function fill_workspace_arrays!(line, z_convs, grid, tloop, data_inds, rot_shifts,
                                λΔDs, wavall, widall, lwavgrid, rwavgrid)
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
            # find position on disk and move to next iter if off disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
            if r2 > 1.0
                @inbounds λΔDs[i,j] = 0.0
                continue
            end

            # calculate the shifted center of the line
            @inbounds λΔDs[i,j] = line * (1.0 + rot_shifts[i,j]) * (1.0 + z_convs)

            # slice out the correct views of the input data for position
            wavt = CUDA.view(wavall, :, tloop[i,j], data_inds[i,j])
            widt = CUDA.view(widall, :, tloop[i,j], data_inds[i,j])

            len = CUDA.size(rwavgrid,3)
            for k in idz:sdz:CUDA.size(rwavgrid,3)
                # set tgrids based on bisector + wid data
                @inbounds lwavgrid[i,j,k] = (λΔDs[i,j] - (0.5 * widt[k] - wavt[k]))
                @inbounds rwavgrid[i,j,k] = (λΔDs[i,j] + (0.5 * widt[k] + wavt[k]))
                if k == 1
                    @inbounds rwavgrid[i,j,1] = lwavgrid[i,j,1] + 1e-3
                end
            end
        end
    end
    return nothing
end

function concatenate_workspace_arrays!(grid, tloop, data_inds, depall,
                                       lwavgrid, rwavgrid, allwavs, allints)
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
            # find position on disk and move to next iter if off disk
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
            if r2 > 1.0
                continue
            end

            # slice out the correct views of the input data for position
            dept = CUDA.view(depall, :, tloop[i,j], data_inds[i,j])

            len = CUDA.size(rwavgrid,3)
            for k in idz:sdz:CUDA.size(rwavgrid,3)
                @inbounds allwavs[i,j,k+len] = rwavgrid[i,j,k]
                @inbounds allints[i,j,k+len] = dept[k]
                @inbounds allwavs[i,j,k] = lwavgrid[i,j, len - (k - 1)]
                @inbounds allints[i,j,k] = dept[len - (k - 1)]
            end
        end
    end
    return nothing
end


function line_profile_gpu!(star_map, grid, lambdas, λΔDs, allwavs, allints)
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
            # set intensity to zero and go to next iter if off grid
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
            if r2 > 1.0
                continue
            end

            # take view of arrays to pass to interpolater
            allwavs_ij = CUDA.view(allwavs, i, j, :)
            allints_ij = CUDA.view(allints, i, j, :)

            for k in idz:sdz:CUDA.length(lambdas)
                # skip all this work if far from line core
                # if (lambdas[k] < (λΔDs[i,j] - 1.0)) || (lambdas[k] > (λΔDs[i,j] + 1.0))
                #     @cuprintln("derp")
                #     continue
                # end

                # do the interpolation
                if ((lambdas[k] < CUDA.first(allwavs_ij)) || (lambdas[k] > CUDA.last(allwavs_ij)))
                    # factor = 1.0
                    continue
                else
                    m = CUDA.searchsortedfirst(allwavs_ij, lambdas[k]) - 1
                    m0 = CUDA.clamp(m, CUDA.firstindex(allints_ij), CUDA.lastindex(allints_ij))
                    m1 = CUDA.clamp(m+1, CUDA.firstindex(allints_ij), CUDA.lastindex(allints_ij))
                    # factor = (allints_ij[m0] * (allwavs_ij[m1] - lambdas[k]) + allints_ij[m1] * (lambdas[k] - allwavs_ij[m0])) / (allwavs_ij[m1] - allwavs_ij[m0])
                    @inbounds star_map[i,j,k] *= ((allints_ij[m0] * (allwavs_ij[m1] - lambdas[k]) + allints_ij[m1] * (lambdas[k] - allwavs_ij[m0])) / (allwavs_ij[m1] - allwavs_ij[m0]))
                end
                # @inbounds star_map[i,j,k] = star_map[i,j,k] * factor
            end
        end
    end
    return nothing
end
