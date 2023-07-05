function fill_workspaces!(line, variability, extra_z, grid, tloop, data_inds, z_rot,
                          z_cbs, bisall, intall, widall, allwavs, allints)
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

            # calculate shifted line center
            λΔD = line * (1.0 + z_rot[i,j]) * (1.0 + z_cbs[i,j] * variability) * (1.0 + extra_z)

            # get length of input data arrays to loop over
            lent = 100
            for k in idz:sdz:lent
                # get forward and reverse indices
                idx1 = k
                idx2 = lent - (k - 1)

                # slice out the correct views of the input data for position
                @inbounds bis1 = bisall[idx1, tloop[i,j], data_inds[i,j]]
                @inbounds wid1 = widall[idx1, tloop[i,j], data_inds[i,j]]
                @inbounds int1 = intall[idx1, tloop[i,j], data_inds[i,j]]

                @inbounds bis2 = bisall[idx2, tloop[i,j], data_inds[i,j]]
                @inbounds wid2 = widall[idx2, tloop[i,j], data_inds[i,j]]
                @inbounds int2 = intall[idx2, tloop[i,j], data_inds[i,j]]

                # right side of line, indexing from middle left to right
                @inbounds allwavs[i,j,k+lent] = (λΔD + (0.5 * wid1 + bis1))
                @inbounds allints[i,j,k+lent] = int1

                # left sight of line, indexing from middle right to left
                @inbounds allwavs[i,j,k] = (λΔD - (0.5 * wid2 - bis2))
                @inbounds allints[i,j,k] = int2
            end
        end
    end
    return nothing
end


function line_profile_gpu!(star_map, grid, lambdas, allwavs, allints)
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

            # set up interpolator
            itp = linear_interp_gpu(allwavs_ij, allints_ij)

            # loop over wavelengths
            for k in idz:sdz:CUDA.length(lambdas)
                if ((lambdas[k] < CUDA.first(allwavs_ij)) || (lambdas[k] > CUDA.last(allwavs_ij)))
                    continue
                else
                    @inbounds star_map[i,j,k] *= itp(lambdas[k])
                end
            end
        end
    end
    return nothing
end
