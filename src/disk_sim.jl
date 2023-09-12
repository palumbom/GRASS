function disk_sim(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                  wsp::SynthWorkspace{T}, prof::AA{T,1}, flux::AA{T,2},
                  tloop::AA{Int,1}; verbose::Bool=true,
                  skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # get sum of weights
    sum_wts = sum(wsp.wts)
    z_cbs_avg = sum(wsp.wts .* wsp.cbs) / sum_wts

    # loop over time
    for t in 1:disk.Nt
        # if skip times is true, continue to next iter
        if skip_times[t]
            tloop .+= 1
            continue
        end

        # loop over wavelength
        for l in 1:length(spec.lines)
            # reset prof
            prof .= zero(T)

            # loop over spatial patches
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

                # check that tloop hasn't exceeded number of epochs
                if tloop[i] > len
                    tloop[i] -= len
                end

                # get views needed for line synthesis
                wsp.bist .= copy(view(soldata.bis[key], :, tloop[i]))
                wsp.intt .= copy(view(soldata.int[key], :, tloop[i]))
                wsp.widt .= copy(view(soldata.wid[key], :, tloop[i]))

                # get amount of convective blueshift needed
                extra_z = spec.conv_blueshifts[l] - z_cbs_avg

                # get shifted line center
                λΔD = spec.lines[l]
                λΔD *= (1.0 + z_rot)
                λΔD *= (1.0 + z_cbs .* spec.variability[l])
                λΔD *= (1.0 + extra_z)

                # get rid of bisector and fix width if variability is turned off
                wsp.bist .*= spec.variability[l]
                if !spec.variability[l]
                    wsp.widt .= view(soldata.wid[key], :, 1)
                end

                # get depth to trim to from depth contrast
                dtrim = spec.depths[l] * soldata.dep_contrast[key]

                # first trim the bisectors to the correct depth
                trim_bisector!(dtrim, wsp.bist, wsp.intt)

                # update the line profile in place
                line_profile_cpu!(λΔD, wsp.wts[i], spec.lambdas, prof, wsp)
            end

            # apply normalization term and add to flux
            flux[:,t] .*= prof./ sum_wts
        end

        # iterate tloop
        tloop .+= 1
    end

    # set instances of outspec where skip is true to 0 and return
    flux[:, skip_times] .= zero(T)
    return nothing
end

function calc_rossiter_quantities!(xyz_planet::AA{T,1}, planet::Planet{T},
                                   disk::DiskParams{T}, wsp::SynthWorkspace{T},
                                   ros_allocs::RossiterAllocs{T}) where T<:AF
    # parse out composite type fields
    Nsubgrid = disk.Nsubgrid

    # allocate memory that wont be needed outside this function
    d2_sub = zeros(Nsubgrid, Nsubgrid)
    μs_sub = zeros(Nsubgrid, Nsubgrid)
    ld_sub = zeros(Nsubgrid, Nsubgrid)
    dA_sub = zeros(Nsubgrid, Nsubgrid)
    dp_sub = zeros(Nsubgrid, Nsubgrid)
    xyz_sub = repeat([zeros(3)], Nsubgrid, Nsubgrid)
    z_rot_sub = zeros(Nsubgrid, Nsubgrid)
    idx = BitMatrix(undef, size(μs_sub))
    idx1 = BitMatrix(undef, size(μs_sub))
    idx2 = BitMatrix(undef, size(μs_sub))

    # loop over disk positions
    for i in eachindex(wsp.ϕc)
        # get the number of theta tiles needed for the latitude tiles
        Nθ = get_Nθ(wsp.ϕc[i], step(disk.ϕe))

        # get the latitude and longitude increments
        dϕ = step(disk.ϕe)
        dθ = deg2rad(360.0) / Nθ

        # get edges of large tile
        ϕ_l = wsp.ϕc[i] - dϕ/2.0
        ϕ_r = wsp.ϕc[i] + dϕ/2.0
        θ_l = wsp.θc[i] + dθ/2.0
        θ_r = wsp.θc[i] + dθ/2.0

        # subdivide the tile
        ϕe_sub = range(ϕ_l, ϕ_r, length=Nsubgrid+1)
        θe_sub = range(θ_l, θ_r, length=Nsubgrid+1)
        ϕc_sub = get_grid_centers(ϕe_sub)
        θc_sub = get_grid_centers(θe_sub)
        subgrid = Iterators.product(ϕc_sub, θc_sub)

        # get cartesian coord for each subgrid and rotate by rot. matrix
        xyz_sub .= map(x -> sphere_to_cart.(disk.ρs, x...), subgrid)
        xyz_sub .= map(x -> disk.R_x * x, xyz_sub)

        # calculate mu at each point
        μs_sub .= map(x -> calc_mu(x, disk.O⃗), xyz_sub)

        # calculate the distance between subtile center and planet
        d2_sub .= map(x -> calc_proj_dist2(x, xyz_planet), xyz_sub)

        # if entire course tile visible, use old weights and move on
        if all(d2_sub .> planet.radius^2.0)
            continue
        end

        # calc limb darkening
        ld_sub .= map(x -> quad_limb_darkening(x, disk.u1, disk.u2), μs_sub)

        # get rotational velocity for location on disk
        z_rot_sub .= map(x -> patch_velocity_los(x..., disk), subgrid)

        # calculate area element of tile
        dA_sub .= map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕe_sub), step(θe_sub)), subgrid)
        dp_sub .= map(x -> abs(dot(x .- disk.O⃗, x)), xyz_sub) / norm(disk.O⃗)

        # get indices for visible patches
        idx1 = μs_sub .> 0.0
        idx2 = d2_sub .> planet.radius^2.0
        idx .= (idx1 .& idx2)

        # if no patches are zero, set wts, etc. to zero and move on
        if all(iszero(idx))
            ros_allocs.μs[i] = 0.0
            ros_allocs.ld[i] = 0.0
            ros_allocs.dA[i] = 0.0
            ros_allocs.wts[i] = 0.0
            ros_allocs.z_rot[i] = 0.0
            continue
        end

        # get total projected, visible area of larger tile
        dA_total = sum(view(dA_sub, idx))
        dA_total_proj = sum(view(dA_sub .* dp_sub, idx))

        # set limb darkening as mean of visible patches
        ros_allocs.μs[i] = mean(view(μs_sub, idx))
        ros_allocs.ld[i] = mean(view(ld_sub, idx))
        ros_allocs.dA[i] = dA_total_proj

        ros_allocs.wts[i] = mean(view(ld_sub .* dA_total_proj, idx)) / wsp.wts[i]
        ros_allocs.z_rot[i] = mean(view(z_rot_sub, idx))
    end

    return nothing
end

function disk_sim_rossiter(spec::SpecParams{T}, disk::DiskParams{T}, planet::Planet{T},
                           soldata::SolarData{T}, wsp::SynthWorkspace{T}, ros_allocs::RossiterAllocs{T},
                           prof::AA{T,1}, flux::AA{T,2}, tloop::AA{Int,1}; verbose::Bool=true,
                           skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # get sum of weights
    sum_wts = sum(wsp.wts)
    z_cbs_avg = sum(wsp.wts .* wsp.cbs) / sum_wts

    # loop over time
    for t in 1:disk.Nt
        # if skip times is true, continue to next iter
        if skip_times[t]
            tloop .+= 1
            continue
        end

        # recopy weights, etc. from unobstructed disk
        ros_allocs.μs .= wsp.μs
        ros_allocs.ld .= wsp.ld
        ros_allocs.dA .= wsp.dA
        ros_allocs.wts .= wsp.wts
        ros_allocs.z_rot .= wsp.z_rot

        # get the projected position of the center of the planet
        # xyz_planet = calc_projected_planet_position()
        xyz_planet = [25.0, 0.0, 1.25]
        xyz_star = [0.0, 0.0, 0.0]

        if t == 2
            xyz_planet = [0.5, 0.0, 1.25]
        end

        # calculate projected distance between planet and star centers
        dist2 = calc_proj_dist2(xyz_planet, xyz_star)

        # check if the planet is transiting
        if dist2 <= (disk.ρs + planet.radius)^2.0
            # re-calculate patch weights, etc. for occulted patches
            calc_rossiter_quantities!(xyz_planet, planet, disk, wsp, ros_allocs)
        end

        # loop over wavelength
        for l in 1:length(spec.lines)
            # reset prof
            prof .= zero(T)

            # loop over spatial patches
            for i in eachindex(wsp.μs)
                # move to next patch if entire patch is behind hemisphere
                μc = ros_allocs.μs[i]
                μc <= zero(T) && continue

                # move to next iteration if entire patch is occulted
                wts = ros_allocs.wts[i] * wsp.wts[i]
                # iszero(wts) && continue

                # get input data for place on disk
                key = wsp.keys[i]
                len = soldata.len[key]

                # get total desired convective blueshift for line
                z_cbs = wsp.cbs[i]

                # get rotational shift
                z_rot = ros_allocs.z_rot[i]

                # check that tloop hasn't exceeded number of epochs
                if tloop[i] > len
                    tloop[i] -= len
                end

                # get views needed for line synthesis
                wsp.bist .= copy(view(soldata.bis[key], :, tloop[i]))
                wsp.intt .= copy(view(soldata.int[key], :, tloop[i]))
                wsp.widt .= copy(view(soldata.wid[key], :, tloop[i]))

                # get amount of convective blueshift needed
                extra_z = spec.conv_blueshifts[l] - z_cbs_avg
                @show extra_z

                # get shifted line center
                λΔD = spec.lines[l]
                λΔD *= (1.0 + z_rot)
                λΔD *= (1.0 + z_cbs .* spec.variability[l])
                λΔD *= (1.0 + extra_z)

                # get rid of bisector and fix width if variability is turned off
                wsp.bist .*= spec.variability[l]
                if !spec.variability[l]
                    wsp.widt .= view(soldata.wid[key], :, 1)
                end

                # get depth to trim to from depth contrast
                dtrim = spec.depths[l] * soldata.dep_contrast[key]

                # first trim the bisectors to the correct depth
                trim_bisector!(dtrim, wsp.bist, wsp.intt)

                # update the line profile in place
                line_profile_cpu!(λΔD, wts, spec.lambdas, prof, wsp)
            end

            # apply normalization term and add to flux
            flux[:,t] .*= prof ./ sum_wts
        end

        # iterate tloop
        tloop .+= 1
    end

    # set instances of outspec where skip is true to 0 and return
    flux[:, skip_times] .= zero(T)
    return nothing
end
