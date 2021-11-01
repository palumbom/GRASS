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

function disk_sim(star_map, tstart, lines, depths, z_convs, grid, lambdas,
                  disc_ax, disc_mu, lenall, wavall, bisall, widall, depall,
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
            # set intensity to zero and go to next iter if off grid
            x = grid[i]
            y = grid[j]
            r2 = calc_r2(x, y)
            if r2 > 1.0
                star_map[i,j,:] .= 0.0
                continue
            end

            # calculate mu for limb darkening
            mu = calc_mu(r2)

            # get rotational redshift
            z_rot = patch_velocity_los(x, y)

            # find the nearest mu ind and ax code
            nn_mu_ind = searchsortednearest_gpu(disc_mu, mu)
            nn_ax_code = find_nearest_ax_gpu(x, y)

            # find the correct data index
            # data_ind = find_data_index_gpu(nn_mu_ind, nn_ax_code)
            data_ind = 41

            # find out number of time epochs of input data for position
            ntimes = lenall[data_ind]
            tloop = tstart[i,j]
            if tloop > ntimes
                tloop = 1
                @inbounds tstart[i,j] = 1
            end

            # slice out the correct views
            wavt = CUDA.view(wavall, :, tloop, data_ind)
            bist = CUDA.view(bisall, :, tloop, data_ind)
            widt = CUDA.view(widall, :, tloop, data_ind)
            dept = CUDA.view(depall, :, tloop, data_ind)

            # iterate tstart
            @inbounds tstart[i,j] += 1

            # loop over lines
            for k in idz:sdz:CUDA.length(lambdas)
                # calculate limb darkening
                @inbounds star_map[i,j,k] = calc_norm_term(mu, CUDA.length(grid), 0.4, 0.26)

                # loop over lines
                for l in 1:CUDA.length(lines)
                    # trim the input data
                    trim_bisector_chop_gpu!(depths[l], wavt, bist, dept, widt, NaN)

                    # calculate the shifted center of the line
                    λΔD = lines[l] * (1.0 + z_rot) * (1.0 + z_convs[l])

                    # skip all this work if far from line core
                    # if (lambdas[k] < (λΔD - 0.5)) | (lambdas[k] > (λΔD + 0.5))
                    #     continue
                    # end

                    # set wavgrids to line center to start
                    for n in 1:CUDA.size(lwavgrid,3)
                        @inbounds lwavgrid[i,j,n] = (λΔD - (0.5 * widt[n] - wavt[n]))
                        @inbounds rwavgrid[i,j,n] = (λΔD + (0.5 * widt[n] + wavt[n]))
                    end
                    @inbounds rwavgrid[i,j,1] = lwavgrid[i,j,1] + 1e-3

                    # concatenate into one big array
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

                    # interpolate onto original lambda grid, extrapolate to continuum
                    factor = linear_interp_mult_gpu(lambdas[k], allwavs_ij, allints_ij, 1.0)
                    @inbounds star_map[i,j,k] = factor * star_map[i,j,k]
                end
            end
        end
    end
    return nothing
end
