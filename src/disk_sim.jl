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

function precompute_quantities_old(wsp::SynthWorkspace{T}, disk::DiskParams{T}, soldata::SolarData{T}) where T<:AF
    # calculate normalization terms and get convective blueshifts
    numer = 0
    denom = 0
    xyz = zeros(3)

    # get discrete mu and ax values
    disc_mu = soldata.mu
    disc_ax = soldata.ax

    # loop over disk positions
    for i in eachindex(disk.ϕc)
        for j in eachindex(disk.θc)
            # get cartesian coord
            xyz .= sphere_to_cart(disk.ρs, disk.ϕc[i], disk.θc[j])
            rotate_vector!(xyz, disk.R_θ)

            # calculate mu
            μc = calc_mu(xyz, disk.O⃗)
            wsp.μs[i,j] = μc

            # move to next iteration if patch element is not visible
            μc <= 0.0 && continue

            # get input data for place on disk
            key = get_key_for_pos(μc, xyz[1], xyz[3], disc_mu, disc_ax)

            # calc limb darkening normalization term
            ld = quad_limb_darkening(μc, disk.u1, disk.u2)
            dA = calc_projected_dA(disk.ϕc[i], disk.θc[j], disk)
            numer += soldata.cbs[key] * (ld * dA)
            denom += ld * dA

            # get rotational velocity for location on disk
            z_rot = patch_velocity_los(disk.ϕc[i], disk.θc[j], disk)

            # copy to workspace
            wsp.dA[i,j] = dA
            wsp.ld[i,j] = ld
            wsp.z_rot[i,j] = z_rot
            wsp.keys[i,j] = key
        end
    end
    return numer/denom, denom
end

function precompute_quantities(xyz::Matrix{Vector{T}}, wsp::SynthWorkspace{T}, disk::DiskParams{T},
                               soldata::SolarData{T}; Nsubgrid::Int=50) where T<:AF
    # calculate normalization terms and get convective blueshifts
    numer = 0
    denom = 0
    # xyz = zeros(Nsubgrid, Nsubgrid, 3)

    # get discrete mu and ax values
    disc_mu = soldata.mu
    disc_ax = soldata.ax

    # loop over disk positions
    for i in eachindex(disk.ϕc)
        for j in 1:disk.Nθ[i]
            # subdivide the tile
            ϕsub = range(disk.ϕe[i], disk.ϕe[i+1], length=Nsubgrid)
            θsub = range(disk.θe[i,j], disk.θe[i,j+1], length=Nsubgrid)
            subgrid = Iterators.product(ϕsub, θsub)

            # get cartesian coord for each subgrid and rotate by rot. matrix
            xyz .= map(x -> sphere_to_cart.(disk.ρs, x...), subgrid)
            xyz .= map(x -> disk.R_θ * x, xyz)

            # calculate mu at each point
            μs = map(x -> calc_mu(x, disk.O⃗), xyz)

            # move to next iteration if patch element is not visible
            all(μs .<= zero(T)) && continue

            # assign the mean mu as the mean of visible mus
            idx = μs .> 0.0
            wsp.μs[i,j] = mean(view(μs, idx))

            # find xyz at mean value of mu
            mean_x = mean(view(getindex.(xyz,1), idx))
            mean_z = mean(view(getindex.(xyz,3), idx))

            # get input data for place on disk
            key = get_key_for_pos(wsp.μs[i,j], mean_x, mean_z, disc_mu, disc_ax)

            # calc limb darkening
            ld = map(x -> quad_limb_darkening(x, disk.u1, disk.u2), μs)

            # get rotational velocity for location on disk
            z_rot = map(x -> patch_velocity_los(x..., disk), subgrid)

            # calculate area element of tile
            dA = map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕsub), step(θsub)), subgrid)
            dA .*= map(x -> abs(dot(x .- disk.O⃗, x)), xyz)

            # copy to workspace
            wsp.dA[i,j] = mean(dA[idx])
            wsp.ld[i,j] = mean(ld[idx])
            wsp.z_rot[i,j] = mean(z_rot[idx])
            wsp.keys[i,j] = key

            # get disk-averaged cbs
            numer += soldata.cbs[key] * (wsp.ld[i,j] * wsp.dA[i,j])
            denom += wsp.ld[i,j] * wsp.dA[i,j]
        end
    end
    return numer/denom, denom
end

function disk_sim_3d(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                     wsp::SynthWorkspace, prof::AA{T,1}, outspec::AA{T,2},
                     tloop::AA{Int,2}; verbose::Bool=true,
                     skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # set pre-allocations and make generator that will be re-used
    outspec .= zero(T)
    liter = 1:length(spec.lines); @assert length(liter) >= 1

    # get intensity-weighted disk-avereged convective blueshift
    z_cbs_avg, sum_norm_terms = precompute_quantities(wsp, disk, soldata)

    # loop over grid positions
    for i in eachindex(disk.ϕc)
        for j in eachindex(disk.θc)
            # move to next iteration if patch element is not visible
            μc = wsp.μs[i,j]
            μc <= zero(T) && continue

            # get input data for place on disk
            key = wsp.keys[i,j]
            len = soldata.len[key]

            # get total desired convective blueshift for line
            z_cbs = soldata.cbs[key]

            # get ld and projected area element
            ld = wsp.ld[i,j]
            dA = wsp.dA[i,j]
            z_rot = wsp.z_rot[i,j]

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
                time_loop_cpu(tloop[i,j], prof, z_rot, z_cbs, z_cbs_avg,
                              key, liter, spec, soldata, wsp)

                # apply normalization term and add to outspec
                outspec[:,t] .+= (prof .* ld * dA)

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
