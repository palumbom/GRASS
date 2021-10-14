"""
    synthesize_spectra(spec, disk; seed_rng=false, verbose=true, top=NaN)

Synthesize spectra given parameters in `spec` and `disk` instances.

# Arguments
- `spec::SpecParams`: SpecParams instance
- `disk::DiskParams`: DiskParams instance
"""
function synthesize_spectra(spec::SpecParams, disk::DiskParams;
                            seed_rng::Bool=false, verbose::Bool=true,
                            top::Float64=NaN)
    # figure out dims for allocation
    if use_gpu
        dims = [length(spec.lines)]
    else
        dims = []
    end

    # allocate memory for synthesis
    Nλ = length(spec.lambdas)
    # prof = ArrayType(ones(Nλ, dims...))
    prof = ArrayType(ones(Nλ))
    outspec = zeros(Nλ, disk.Nt)

    # run the simulation (outspec modified in place)
    disk_sim(spec, disk, prof, outspec, seed_rng=seed_rng, verbose=verbose, top=top)
    return Array(spec.lambdas), outspec
end

# line loop function, update prof in place
function line_loop_cpu(prof::AA{T,1}, mid::T, depth::T, rot_shift::T,
                       conv_blueshift::T, lambdas::AA{T,1},
                       wsp::SynthWorkspace{T}; top::T=NaN) where T<:AF
    # first trim the bisectors to the correct depth
    trim_bisector_chop!(depth, wsp.wavt, wsp.bist, wsp.dept, wsp.widt, top=top)

    # calculate line center given rot. and conv. doppler shift -> λrest * (1 + z)
    λΔD = mid * (one(T) + rot_shift) * (one(T) + conv_blueshift)

    # find window around shifted line
    lind = findfirst(x -> x > λΔD - 0.5, lambdas)
    rind = findfirst(x -> x > λΔD + 0.5, lambdas)
    lambda_window = view(lambdas, lind:rind)
    prof_window = view(prof, lind:rind)

    # update the line profile in place
    line_profile!(λΔD, lambda_window, prof_window, wsp)
    return nothing
end

function time_loop_cpu(t_loop::Int, prof::AA{T,1}, rot_shift::T,
                       key::Tuple{Symbol, Symbol}, liter::UnitRange{Int},
                       spec::SpecParams{T}, wsp::SynthWorkspace{T}; top::T=NaN) where T<:AF
    # some assertions
    @assert all(prof .== one(T))

    # get views needed for line synthesis
    wsp.wavt .= view(spec.soldata.wav[key], :, t_loop)
    wsp.bist .= view(spec.soldata.bis[key], :, t_loop)
    wsp.dept .= view(spec.soldata.dep[key], :, t_loop)
    wsp.widt .= view(spec.soldata.wid[key], :, t_loop)

    # loop over specified synthetic lines
    for l in liter
        wsp.wavt .*= spec.variability[l]
        line_loop(prof, spec.lines[l], spec.depths[l], rot_shift,
                  spec.conv_blueshifts[l], spec.lambdas, wsp, top=top)
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

function disk_sim(spec::SpecParams{T}, disk::DiskParams{T,Int64}, prof::AA{T,1},
                  outspec::AA{T,2}; top::T=NaN, seed_rng::Bool=false,
                  skip_times::BitVector=BitVector(zeros(disk.Nt)),
                  verbose::Bool=true) where T<:AF
    # make grid
    grid = make_grid(N=disk.N)

    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    wsp = SynthWorkspace(spec)
    liter = 1:length(spec.lines)

    # seeding rng
    if seed_rng
        if verbose println("Seeding RNG") end
        Random.seed!(42)
    end

    # loop over grid positions
    for i in grid
        for j in grid
            # move to next iteration if off grid
            (i^2 + j^2) > one(T) && continue

            # get input data for place on disk
            key = get_key_for_pos(i, j)
            len = spec.soldata.len[key]

            # get redshift z and norm term for location on disk
            rot_shift = patch_velocity_los(i, j, pole=disk.pole)
            norm_term = calc_norm_term(i, j, disk)

            # loop over time, starting at random epoch
            inds = generate_indices(disk.Nt, len)
            for (t, t_loop) in enumerate(inds)
                # if skip times is true, continue to next iter
                skip_times[t] && continue

                # update profile in place
                prof .= one(T)
                time_loop(t_loop, prof, rot_shift, key, liter, spec, wsp, top=top)

                # apply normalization term and add to outspec
                outspec[:,t] .+= (Array(prof) .* norm_term)
            end
        end
    end

    # set instances of outspec where skip is true to 0 and return
    outspec ./= maximum(outspec, dims=1)
    outspec[:, skip_times] .= zero(T)
    return nothing
end
