# line loop function, update prof in place
function line_loop_cpu(prof::AA{T,1}, λΔD::T, depth::T, lambdas::AA{T,1},
                       wsp::SynthWorkspace{T}) where T<:AF
    # first trim the bisectors to the correct depth
    trim_bisector!(depth, wsp.bist, wsp.intt)

    # update the line profile in place
    line_profile_cpu!(λΔD, lambdas, prof, wsp)
    return nothing
end

function time_loop_cpu(tloop::Int, prof::AA{T,1}, z_rot::T, z_cbs::T,
                       z_cbs_avg::T, key::Tuple{Symbol, Symbol},
                       liter::UnitRange{Int}, spec::SpecParams{T},
                       soldata::SolarData, wsp::SynthWorkspace{T}) where T<:AF
    # reset prof
    prof .= one(T)

    # loop over lines
    for l in liter
        # get views needed for line synthesis
        wsp.bist .= copy(view(soldata.bis[key], :, tloop))
        wsp.intt .= copy(view(soldata.int[key], :, tloop))
        wsp.widt .= copy(view(soldata.wid[key], :, tloop))

        # calculate the position of the line center
        extra_z = spec.conv_blueshifts[l] - z_cbs_avg
        λΔD = spec.lines[l] * (1.0 + z_rot) * (1.0 + z_cbs) * (1.0 + extra_z)

        # get rid of bisector and fix width if variability is turned off
        wsp.bist .*= spec.variability[l]
        if !spec.variability[l]
            wsp.widt .= view(soldata.wid[key], :, 1)
        end

        # get depth to trim to from depth contrast
        dtrim = spec.depths[l] * soldata.dep_contrast[key]

        # synthesize the line
        line_loop_cpu(prof, λΔD, dtrim, spec.lambdas, wsp)
    end
    return nothing
end

function calc_disk_avg_cbs(disk::DiskParams{T}, soldata::SolarData{T},
                           grid::StepRangeLen, disc_mu::AA{T,1},
                           disc_ax::AA{Int,1}) where T<:AF
    # calculate normalization terms and get convective blueshifts
    numer = 0
    denom = 0
    for i in eachindex(grid)
        for j in eachindex(grid)
            # get positiosns
            x = grid[i]
            y = grid[j]

            # move to next iteration if off grid
            (x^2 + y^2) > one(T) && continue

            # get input data for place on disk
            key = get_key_for_pos(x, y, disc_mu, disc_ax)

            # calc limb darkening and get convective blueshift
            norm_term = calc_norm_term(x, y, disk)
            numer += soldata.cbs[key] * norm_term
            denom += norm_term
        end
    end
    return numer/denom, denom
end

function disk_sim(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                  prof::AA{T,1}, outspec::AA{T,2}, tloop::AA{Int,2}; verbose::Bool=true,
                  skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # make grid
    grid = make_grid(N=disk.N)

    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    wsp = SynthWorkspace()
    liter = 1:length(spec.lines); @assert length(liter) >= 1

    # get the value of mu and ax codes
    disc_ax = parse_ax_string.(getindex.(keys(soldata.len),1))
    disc_mu = parse_mu_string.(getindex.(keys(soldata.len),2))

    # get indices to sort by mus
    inds_mu = sortperm(disc_mu)
    disc_mu .= disc_mu[inds_mu]
    disc_ax .= disc_ax[inds_mu]

    # get indices to sort by axis within mu sort
    for mu_val in unique(disc_mu)
        inds1 = (disc_mu .== mu_val)
        inds2 = sortperm(disc_ax[inds1])

        disc_mu[inds1] .= disc_mu[inds1][inds2]
        disc_ax[inds1] .= disc_ax[inds1][inds2]
    end

    # get intensity-weighted disk-avereged convective blueshift
    z_cbs_avg, sum_norm_terms = calc_disk_avg_cbs(disk, soldata, grid, disc_mu, disc_ax)

    # loop over grid positions
    for i in eachindex(grid)
        for j in eachindex(grid)
            # get positiosns
            x = grid[i]
            y = grid[j]

            # move to next iteration if off grid
            (x^2 + y^2) > one(T) && continue

            # get input data for place on disk
            key = get_key_for_pos(x, y, disc_mu, disc_ax)
            len = soldata.len[key]

            # get total doppler shift for the line, and norm_term
            z_cbs = soldata.cbs[key]
            z_rot = patch_velocity_los(x, y, pole=disk.pole)
            norm_term = calc_norm_term(x, y, disk)

            # loop over time
            for t in 1:disk.Nt
                # check that tloop hasn't exceeded number of epochs
                if tloop[i,j] > len
                    tloop[i,j] = 1
                end

                # if skip times is true, continue to next iter
                if skip_times[t]
                    tloop[i,j] += 1
                    continue
                end

                # update profile in place
                time_loop_cpu(tloop[i,j], prof, z_rot, z_cbs, z_cbs_avg, key, liter, spec, soldata, wsp)

                # apply normalization term and add to outspec
                outspec[:,t] .+= (prof .* norm_term)

                # iterate tloop
                tloop[i,j] += 1
            end
        end
    end

    # divide by sum of weights
    outspec ./= sum_norm_terms

    # set instances of outspec where skip is true to 0 and return
    outspec[:, skip_times] .= zero(T)
    return nothing
end
