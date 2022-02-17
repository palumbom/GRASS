# line loop function, update prof in place
function line_loop_cpu(prof::AA{T,1}, mid::T, depth::T, z_rot::T,
                       conv_blueshift::T, lambdas::AA{T,1},
                       wsp::SynthWorkspace{T}) where T<:AF
    # first trim the bisectors to the correct depth
    trim_bisector!(depth, wsp.wavt, wsp.bist, wsp.dept, wsp.widt)

    # calculate line center given rot. and conv. doppler shift -> λrest * (1 + z)
    λΔD = mid * (one(T) + z_rot) * (one(T) + conv_blueshift)

    # find window around shifted line
    lind = findfirst(x -> x > λΔD - 1.0, lambdas)
    if isnothing(lind)
        lind = firstindex(lambdas)
    end

    rind = findfirst(x -> x > λΔD + 1.0, lambdas)
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

function time_loop_cpu(t_loop::Int, prof::AA{T,1}, z_rot::T,
                       key::Tuple{Symbol, Symbol}, liter::UnitRange{Int},
                       spec::SpecParams{T}, soldata::SolarData,
                       wsp::SynthWorkspace{T}) where T<:AF
    # get views needed for line synthesis
    wsp.wavt .= view(soldata.wav[key], :, t_loop)
    wsp.bist .= view(soldata.bis[key], :, t_loop)
    wsp.dept .= view(soldata.dep[key], :, t_loop)
    wsp.widt .= view(soldata.wid[key], :, t_loop)

    # loop over specified synthetic lines
    prof .= one(T)
    for l in liter
        wsp.wavt .*= spec.variability[l]
        line_loop_cpu(prof, spec.lines[l], spec.depths[l], z_rot,
                      spec.conv_blueshifts[l], spec.lambdas, wsp)
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

function disk_sim(spec::SpecParams{T}, disk::DiskParams{T,Int64},
                  soldata::SolarData{T}, prof::AA{T,1}, outspec::AA{T,2};
                  seed_rng::Bool=false, verbose::Bool=true, kwargs...) where T<:AF
    # make grid
    grid = make_grid(N=disk.N)

    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    wsp = SynthWorkspace(spec)
    liter = 1:length(spec.lines)

    # get list of discrete mu's in input data
    mu_symb = soldata.mu
    disc_mu = parse_mu_string.(mu_symb)

    # seeding rng
    if seed_rng
        if verbose println("Seeding RNG") end
        Random.seed!(42)
    end

    # loop over spatial grid positions
    rm = false
    if !rm
        f = (t) -> spatial_loop(t[1], t[2], spec, disk, soldata, wsp, prof,
                                outspec, liter, mu_symb, disc_mu; kwargs...)
        map(f, grid)
    else
        grid1D = make_grid_range(disk.N)
        f = (t) -> spatial_loop_rm(t[1], t[2], grid1D, spec, disk,
                                   soldata, wsp, prof, outspec, liter,
                                   mu_symb, disc_mu; kwargs...)
        map(f, CartesianIndices(collect(grid)))
    end

    # ensure normalization and return
    outspec ./= maximum(outspec, dims=1)
    return nothing
end

function spatial_loop(x::T, y::T, spec::SpecParams{T}, disk::DiskParams{T,Int64},
                      soldata::SolarData{T}, wsp::SynthWorkspace{T},
                      prof::AA{T,1}, outspec::AA{T,2}, liter::UnitRange{Int64},
                      mu_symb::AA{Symbol,1}, disc_mu::AA{T,1};
                      skip_times::BitVector=BitVector(zeros(disk.Nt))) where T<:AF
    # move to next iteration if off grid
    calc_r2(x,y) > one(T) && return nothing

    # get input data for place on disk
    key = get_key_for_pos(x, y, disc_mu, mu_symb)

    # use data for same mu from different axis if axis is missing
    while !(key in keys(soldata.len))
        idx = findfirst(key[1] .== soldata.ax)
        if idx == length(soldata.ax)
            idx = 1
        end
        key = (soldata.ax[idx+1], key[2])
    end
    len = soldata.len[key]

    # get redshift z and norm term for location on disk
    z_rot = patch_velocity_los(x, y, pole=disk.pole)
    norm_term = calc_norm_term(x, y, disk)

    # loop over time, starting at random epoch
    inds = generate_indices(disk.Nt, len)
    for (t, t_loop) in enumerate(inds)
        # if skip times is true, continue to next iter
        skip_times[t] && continue

        # update profile in place
        time_loop_cpu(t_loop, prof, z_rot, key, liter, spec, soldata, wsp)

        # apply normalization term and add to outspec
        outspec[:,t] .+= (prof .* norm_term)
    end
    return nothing
end

function spatial_loop_rm(i::Int64, j::Int64, grid::AA{T,1}, spec::SpecParams{T},
                         disk::DiskParams{T,Int64}, soldata::SolarData{T},
                         wsp::SynthWorkspace{T}, prof::AA{T,1},
                         outspec::AA{T,2}, liter::UnitRange{Int64},
                         mu_symb::AA{Symbol,1}, disc_mu::AA{T,1};
                         skip_times::BitVector=BitVector(zeros(disk.Nt))) where T<:AF
    # get positions and move to next iteration if off grid
    x = grid[i]; y = grid[j];
    calc_r2(x,y) > one(T) && return nothing

    # get input data for place on disk
    key = get_key_for_pos(x, y, disc_mu, mu_symb)

    # use data for same mu from different axis if axis is missing
    while !(key in keys(soldata.len))
        idx = findfirst(key[1] .== soldata.ax)
        if idx == length(soldata.ax)
            idx = 1
        end
        key = (soldata.ax[idx+1], key[2])
    end
    len = soldata.len[key]

    # get redshift z and norm term for location on disk
    z_rot = patch_velocity_los(x, y, pole=disk.pole)
    norm_term = calc_norm_term(x, y, disk)

    # loop over time, starting at random epoch
    inds = generate_indices(disk.Nt, len)
    for (t, t_loop) in enumerate(inds)
        # if skip times is true, continue to next iter
        skip_times[t] && continue

        # update profile in place
        time_loop_cpu(t_loop, prof, z_rot, key, liter, spec, soldata, wsp)

        # apply normalization term and add to outspec
        outspec[:,t] .+= (prof .* norm_term)
    end
    return nothing
end
