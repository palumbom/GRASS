function line_profile_gpu!(star_map, tloop, lines, depths, z_convs, grid,
                           lambdas, data_inds, lenall, wavall, bisall, widall,
                           depall, lwavgrid, rwavgrid, allwavs, allints,
                           polex, poley, polez)
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
                @inbounds star_map[i,j,:] .= 0.0
                continue
            end

            # calculate the shifted center of the line
            # TODO: pole implementation
            z_rot = patch_velocity_los_gpu(x, y, 1.0, polex, poley, polez)
            λΔD = lines * (1.0 + z_rot) * (1.0 + z_convs)

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

            for k in idz:sdz:CUDA.length(lambdas)
                # skip all this work if far from line core
                if (lambdas[k] < (λΔD - 0.5)) | (lambdas[k] > (λΔD + 0.5))
                    continue
                end

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
