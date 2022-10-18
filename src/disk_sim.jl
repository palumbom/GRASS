# line loop function, update prof in place
function line_loop_cpu(prof::AA{T,1}, λΔD::T, depth::T, lambdas::AA{T,1},
                       wsp::SynthWorkspace{T}) where T<:AF
    # first trim the bisectors to the correct depth
    trim_bisector!(depth, wsp.bist, wsp.intt)

    # find window around shifted line
    buff = maximum(wsp.widt) / 2.0
    lind = findfirst(x -> x > λΔD - buff, lambdas)
    if isnothing(lind)
        lind = firstindex(lambdas)
    end

    rind = findfirst(x -> x > λΔD + buff, lambdas)
    if isnothing(rind)
        rind = lastindex(lambdas)
    end

    # only compute flux values on window around the shifted line center
    lambda_window = view(lambdas, lind:rind)
    prof_window = view(prof, lind:rind)

    # update the line profile in place
    line_profile_cpu!(λΔD, lambda_window, prof_window, wsp)
    return nothing
end

function time_loop_cpu(tloop::Int, prof::AA{T,1}, z_rot::T, z_cbs::T,
                       z_cbs_avg::T, key::Tuple{Symbol, Symbol},
                       liter::UnitRange{Int}, spec::SpecParams{T},
                       soldata::SolarData, wsp::SynthWorkspace{T}) where T<:AF
    # get views needed for line synthesis
    wsp.bist .= view(soldata.bis[key], :, tloop)
    wsp.intt .= view(soldata.int[key], :, tloop)
    wsp.widt .= view(soldata.wid[key], :, tloop)

    # loop over specified synthetic lines
    prof .= one(T)
    for l in liter
        # calculate the position of the line center
        extra_z = spec.conv_blueshifts[l] - z_cbs_avg
        λΔD = spec.lines[l] * (1.0 + z_rot) * (1.0 + z_cbs) * (1.0 + extra_z)

        # synthesize the line
        wsp.bist .*= spec.variability[l]
        line_loop_cpu(prof, λΔD, spec.depths[l], spec.lambdas, wsp)
    end
    return nothing
end

function calc_disk_avg_cbs(grid::StepRangeLen, disc_mu::AA{T,1}, mu_symb::AA{Symbol,1},
                           disk::DiskParams{T}, soldata::SolarData{T}) where T<:AF
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
            key = get_key_for_pos(x, y, disc_mu, mu_symb)

            # calc limb darkening and get convective blueshift
            norm_term = calc_norm_term(x, y, disk)
            numer += soldata.cbs[key] * norm_term
            denom += norm_term
        end
    end
    return numer/denom, denom
end

function generate_tloop!(tloop::AA{Int,2}, grid::StepRangeLen, disc_mu::AA{T,1},
                         mu_symb::AA{Symbol,1}, soldata::SolarData{T}) where T<:AF
    # make sure dimensions are correct
    @assert size(tloop) == (length(grid), length(grid))

    for i in eachindex(grid)
        for j in eachindex(grid)
            # get positiosns
            x = grid[i]
            y = grid[j]

            # move to next iteration if off grid
            (x^2 + y^2) > one(T) && continue

            # get input data for place on disk
            key = get_key_for_pos(x, y, disc_mu, mu_symb)
            while !(key in keys(soldata.len))
                idx = findfirst(key[1] .== soldata.ax)
                if isnothing(idx) || idx == length(soldata.ax)
                    idx = 1
                end
                key = (soldata.ax[idx+1], key[2])
            end
            len = soldata.len[key]

            tloop[i,j] = floor(Int, rand() * len) + 1
        end
    end
    return nothing
end

function disk_sim(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                  prof::AA{T,1}, outspec::AA{T,2}, tloop::AA{Int,2}; verbose::Bool=true,
                  skip_times::BitVector=BitVector(zeros(disk.Nt))) where T<:AF
    # make grid
    grid = make_grid(N=disk.N)

    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    wsp = SynthWorkspace()
    liter = 1:length(spec.lines); @assert length(liter) >= 1

    # get list of discrete mu's in input data
    mu_symb = soldata.mu
    disc_mu = parse_mu_string.(mu_symb)

    # get intensity-weighted disk-avereged convective blueshift
    z_cbs_avg, sum_norm_terms = calc_disk_avg_cbs(grid, disc_mu, mu_symb, disk, soldata)

    # fill tloop
    # if all(iszero.(tloop))
    generate_tloop!(tloop, grid, disc_mu, mu_symb, soldata)
    # end

    # loop over grid positions
    for i in eachindex(grid)
        for j in eachindex(grid)
            # get positiosns
            x = grid[i]
            y = grid[j]

            # move to next iteration if off grid
            (x^2 + y^2) > one(T) && continue

            # get input data for place on disk
            key = get_key_for_pos(x, y, disc_mu, mu_symb)

            # use data for same mu from different axis if axis is missing
            while !(key in keys(soldata.len))
                idx = findfirst(key[1] .== soldata.ax)
                if isnothing(idx) || idx == length(soldata.ax)
                    idx = 1
                end
                key = (soldata.ax[idx+1], key[2])
            end
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
