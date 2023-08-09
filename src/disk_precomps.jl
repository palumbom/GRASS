function precompute_quantities!(wsp::SynthWorkspace{T}, disk::DiskParams{T}) where T<:AF
    # parse out composite type fields
    Nsubgrid = disk.Nsubgrid

    # allocate memory that wont be needed outside this function
    xyz = repeat([zeros(3)], Nsubgrid, Nsubgrid)

    # loop over disk positions
    for i in eachindex(disk.ϕc)
        for j in 1:disk.Nθ[i]
            # subdivide the tile
            ϕsub = get_grid_centers(range(disk.ϕe[i], disk.ϕe[i+1], length=Nsubgrid+1))
            θsub = get_grid_centers(range(disk.θe[i,j], disk.θe[i,j+1], length=Nsubgrid+1))
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

            # find xz at mean value of mu and get axis code (i.e., N, E, S, W)
            mean_x = mean(view(getindex.(xyz,1), idx))
            mean_z = mean(view(getindex.(xyz,3), idx))
            wsp.ax_codes[i,j] = find_nearest_ax_code(mean_x, mean_z)

            # calc limb darkening
            ld = map(x -> quad_limb_darkening(x, disk.u1, disk.u2), μs)

            # get rotational velocity for location on disk
            z_rot = map(x -> patch_velocity_los(x..., disk), subgrid)

            # calculate area element of tile
            dA = map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕsub), step(θsub)), subgrid)
            dA .*= map(x -> abs(dot(x .- disk.O⃗, x)), xyz)

            # copy to workspace
            wsp.wts[i,j] = mean(view(ld, idx)) * mean(view(dA, idx))
            wsp.z_rot[i,j] = mean(view(z_rot, idx))
        end
    end
    return nothing
end

function generate_tloop!(tloop::AA{Int,2}, wsp::SynthWorkspace{T}, soldata::SolarData{T}) where T<:AF
    generate_tloop!(tloop, wsp.μs, wsp.keys, soldata.len)
    return nothing
end

function generate_tloop!(tloop::AA{Int,2}, μs::AA{T,2}, keys::AA{Tuple{Symbol, Symbol},2}, len::Dict{Tuple{Symbol, Symbol}, Int64}) where T<:AF
    # loop over μ positions
    for i in CartesianIndices(μs)
        # move on if we are off the grid
        μs[i] <= zero(T) && continue

        # get length of input data for place on disk
        maxlen = len[keys[i]]

        # generate random index
        tloop[i] = floor(Int, rand() * maxlen) + 1
    end
    return nothing
end

function get_keys_and_cbs!(wsp::SynthWorkspace{T}, soldata::SolarData{T}) where T<:AF
    get_keys_and_cbs!(wsp.keys, wsp.μs, wsp.cbs, wsp.ax_codes, soldata)
    return nothing
end


function get_keys_and_cbs!(keys::AA{Tuple{Symbol, Symbol},2}, μs::AA{T,2},
                           cbs::AA{T,2}, ax_codes::AA{Int,2}, soldata::SolarData{T}) where T<:AF
    # get the mu and axis codes
    disc_mu = soldata.mu
    disc_ax = soldata.ax

    # loop over μ positions
    for i in CartesianIndices(μs)
        # move on if we are off the grid
        μs[i] <= zero(T) && continue

        # get input data for place on disk
        the_key = get_key_for_pos(μs[i], ax_codes[i], disc_mu, disc_ax)
        cbs[i] = soldata.cbs[the_key]
        keys[i] = the_key
    end
    return nothing
end
