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

function time_loop_cpu(t_loop::Int, prof::AA{T,1}, z_rot::T, z_cbs::T,
                       z_cbs_avg::T, key::Tuple{Symbol, Symbol},
                       liter::UnitRange{Int}, spec::SpecParams{T},
                       soldata::SolarData, wsp::SynthWorkspace{T}) where T<:AF
    # get views needed for line synthesis
    wsp.bist .= view(soldata.bis[key], :, t_loop)
    wsp.intt .= view(soldata.int[key], :, t_loop)
    wsp.widt .= view(soldata.wid[key], :, t_loop)

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

function generate_indices(Nt::Integer, len::Integer)
    # initialize variable to get total number of indices
    sum_lens = 0

    # start at random index less than len and go to len
    start = Iterators.take(rand(1:len):len, Nt)
    sum_lens += length(start)

    # return if that's all that's needed
    if length(start) == Nt
         @assert sum_lens == Nt
        return Iterators.flatten(start)
    end

    # find out how many more cycles are needed
    niter = ceil(Int, (Nt - length(start))/len)

    # make vector of iterators and
    inds = Vector{Base.Iterators.Take{UnitRange{Int64}}}(undef, niter + 1)
    inds[1] = start
    for i in 2:niter
        inds[i] = Iterators.take(1:len, len)
        sum_lens += len
    end

    # ensure the last one only takes the the remainder
    inds[end] = Iterators.take(1:len, Nt - sum_lens)
    sum_lens += length(inds[end])
    @assert sum_lens == Nt

    # return flattened iterator
    return Iterators.flatten(inds)
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

function disk_sim(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                  prof::AA{T,1}, outspec::AA{T,2}; seed_rng::Bool=false,
                  skip_times::BitVector=BitVector(zeros(disk.Nt)),
                  verbose::Bool=true) where T<:AF
    # make grid
    grid = make_grid(N=disk.N)

    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    wsp = SynthWorkspace()
    liter = 1:length(spec.lines); @assert length(liter) >= 1

    # get list of discrete mu's in input data
    mu_symb = soldata.mu
    disc_mu = parse_mu_string.(mu_symb)

    # seeding rng
    if seed_rng
        if verbose println("Seeding RNG") end
        Random.seed!(42)
    end

    # get intensity-weighted disk-avereged convective blueshift
    z_cbs_avg, sum_norm_terms = calc_disk_avg_cbs(grid, disc_mu, mu_symb, disk, soldata)

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

            # get total doppler shift for the line,
            z_cbs = soldata.cbs[key]
            z_rot = patch_velocity_los(x, y, pole=disk.pole)

            # get norm_term for location on disk
            norm_term = calc_norm_term(x, y, disk)

            # loop over time, starting at random epoch
            inds = generate_indices(disk.Nt, len)
            for (t, t_loop) in enumerate(inds)
                # if skip times is true, continue to next iter
                skip_times[t] && continue

                # update profile in place
                time_loop_cpu(t_loop, prof, z_rot, z_cbs, z_cbs_avg, key, liter, spec, soldata, wsp)

                # apply normalization term and add to outspec
                outspec[:,t] .+= (prof .* norm_term)
            end
        end
    end

    # divide by sum of weights
    outspec ./= sum_norm_terms

    # set instances of outspec where skip is true to 0 and return
    outspec[:, skip_times] .= zero(T)
    return nothing
end
