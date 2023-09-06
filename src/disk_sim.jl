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
        λΔD = spec.lines[l] * (1.0 + z_rot) * (1.0 + z_cbs .* spec.variability[l]) * (1.0 + extra_z)

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

function disk_sim(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                  wsp::SynthWorkspace, prof::AA{T,1}, outspec::AA{T,2},
                  tloop::AA{Int,1}; verbose::Bool=true,
                  skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    liter = 1:length(spec.lines); @assert length(liter) >= 1

    # get sum of weights
    sum_wts = sum(wsp.wts)
    z_cbs_avg = sum(wsp.wts .* wsp.cbs) / sum_wts

    # loop over grid position
    for i in eachindex(wsp.μs)
        # move to next iteration if patch element is not visible
        μc = wsp.μs[i]
        μc <= zero(T) && continue

        # get input data for place on disk
        key = wsp.keys[i]
        len = soldata.len[key]

        # get total desired convective blueshift for line
        z_cbs = wsp.cbs[i]

        # get rotational shift
        z_rot = wsp.z_rot[i]

        # loop over time
        for t in 1:disk.Nt
            # check that tloop hasn't exceeded number of epochs
            if tloop[i] > len
                tloop[i] = 1
            end

            # if skip times is true, continue to next iter
            if skip_times[t]
                tloop[i] += 1
                continue
            end

            # update flux profile in place
            time_loop_cpu(tloop[i], prof, z_rot, z_cbs, z_cbs_avg,
                          key, liter, spec, soldata, wsp)

            # apply normalization term and add to outspec
            outspec[:,t] .+= (prof .* wsp.wts[i])

            # iterate tloop
            tloop[i] += 1
        end
    end

    # divide by sum of weights
    outspec ./= sum_wts

    # set instances of outspec where skip is true to 0 and return
    outspec[:, skip_times] .= zero(T)
    return nothing
end

function disk_sim_rm(spec::SpecParams{T}, disk::DiskParams{T}, planet::Planet,
                     soldata::SolarData{T}, prof::AA{T,1}, outspec::AA{T,2},
                     tloop::AA{Int,2}; verbose::Bool=true, nsubgrid::Int=256,
                     skip_times::BitVector=BitVector(zeros(disk.Nt))) where T<:AF

    # pre-calculate the norm terms
    norm_terms = calc_norm_terms(disk)

    # get planet positions
    tvec = get_simulation_times(disk)
    xpos, ypos = calc_planet_position(tvec, planet...)

    # TODO hardcoded for debugging
    xpos = range(-1.25, 1.25, length=disk.Nt)

    # get grid details
    grid_range = make_grid_range(disk.N)
    grid_edges = get_grid_edges(grid_range)

    # allocate memory
    unocculted = trues(nsubgrid, nsubgrid)
    sublimbdarks = zeros(nsubgrid, nsubgrid)
    sub_z_rots = zeros(nsubgrid, nsubgrid)


    # anonymous func for calling spatial loop
    f = (t,n) -> spatial_loop_rm(t[1], t[2], n, grid_range, grid_edges,
                                 spec, disk, planet..., soldata, wsp,
                                 prof, outspec, unocculted, sublimbdarks,
                                 sub_z_rots, liter, mu_symb, disc_mu,
                                 xpos, ypos; nsubgrid=nsubgrid, kwargs...)
    map(f, CartesianIndices((1:disk.N, 1:disk.N)), norm_terms)

    # ensure normalization and return
    outspec ./= maximum(outspec, dims=1)
    return nothing
end

function spatial_loop(x::T, y::T, norm_term::T, spec::SpecParams{T},
                      disk::DiskParams{T}, soldata::SolarData{T},
                      wsp::SynthWorkspace{T}, prof::AA{T,1}, outspec::AA{T,2},
                      liter::UnitRange{Int64}, mu_symb::AA{Symbol,1}, disc_mu::AA{T,1};
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

    # get redshift z for location on disk
    z_rot = patch_velocity_los(x, y, pole=disk.pole)

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

function spatial_loop_rm(i::Int64, j::Int64, norm_term0::T, grid::AA{T,1},
                         grid_edges::AA{T,1}, spec::SpecParams{T},
                         disk::DiskParams{T}, planet::Planet,
                         soldata::SolarData{T}, wsp::SynthWorkspace{T},
                         prof::AA{T,1}, outspec::AA{T,2}, unocculted::AA{Bool,2},
                         sublimbdarks::AA{T,2}, sub_z_rots::AA{T,2},
                         liter::UnitRange{Int64}, mu_symb::AA{Symbol,1},
                         disc_mu::AA{T,1}, xpos::AA{T,1}, ypos::AA{T,1};
                         nsubgrid::Int64=256,
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

    # get initial redshift z for location on disk
    z_rot0 = patch_velocity_los(x, y, pole=disk.pole)

    # loop over time, starting at random epoch
    inds = generate_indices(disk.Nt, len)
    for (t, t_loop) in enumerate(inds)
        # if skip times is true, continue to next iter
        skip_times[t] && continue

        # set z_rot and norm term back to initial value for cell
        z_rot = z_rot0
        norm_term = norm_term0

        # figure out if planet is on patch at this time
        dist2 = calc_dist2(x, y, xpos[t], ypos[t])
        if dist2 <= (planet.radius + step(grid))^2
            xrange = range(grid_edges[i], grid_edges[i+1], length=nsubgrid)
            yrange = range(grid_edges[j], grid_edges[j+1], length=nsubgrid)
            subgrid = Iterators.product(xrange, yrange)

            # calculate subgrid element distances from planet
            f = z -> (calc_dist2(z, (xpos[t], ypos[t])) > planet.radius^2) && calc_mu(z) >= 0.0
            unocculted .= map(f, subgrid)

            # move to next time step if square is completely occluded
            iszero(sum(unocculted)) && continue

            # calculate limb darkening and redshifts in subgrid
            sublimbdarks .= quad_limb_darkening.(calc_mu.(subgrid), disk.u1, disk.u2)
            sub_z_rots .= patch_velocity_los.(subgrid, pole=disk.pole)

            # take views of arrays
            z_rots_unocculted = view(sub_z_rots, unocculted)
            sublimbdarks_unocculted = view(sublimbdarks, unocculted)

            # calculate new weighted z_rot and norm_term
            z_rot = sum(z_rots_unocculted .* sublimbdarks_unocculted) / sum(sublimbdarks_unocculted)
            norm_term *= sum(sublimbdarks_unocculted)/sum(sublimbdarks)
        end

        # update profile in place
        time_loop_cpu(t_loop, prof, z_rot, key, liter, spec, soldata, wsp)

        # apply normalization term and add to outspec
        outspec[:,t] .+= (prof .* norm_term)
    end
    return nothing
end
