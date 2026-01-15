struct SynthWorkspaceEclipse{T<:AF}
    lwavgrid::AA{T,1}
    rwavgrid::AA{T,1}
    allwavs::AA{T,1}
    allints::AA{T,1}

    bist::AA{T,1}
    intt::AA{T,1}
    widt::AA{T,1}

    SP_sun_pos::Matrix{Vector{Float64}} 
    SP_sun_vel::Matrix{Vector{Float64}} 
    v_scalar_grid::AA{T,2} 
    pole_vector_grid::Matrix{Vector{Float64}} 
    SP_bary_pos::Matrix{Vector{Float64}} 
    SP_bary_vel::Matrix{Vector{Float64}} 
    mu_grid::AA{T,2} 
    projected_velocities_no_cb::AA{T,2} 
    SP_bary::Matrix{Vector{Float64}} 
    OP_bary::Matrix{Vector{Float64}} 
    v_earth_orb_proj::Matrix{Float64} 
    distance::AA{T,2} 
    mean_weight_v_no_cb::AA{T,3} 
    mean_weight_v_earth_orb::AA{T,3} 

    z_cbs::AA{T,2}
    ϕc::AA{T,2}
    θc::AA{T,2}
    μs::AA{T,3} 
    ld::AA{T,3} 
    ext::AA{T,3}
    dA::AA{T,3} 
    xyz::AA{T,3} 

    cbs::AA{T,2}
    z_rot::AA{T,3} 
    ax_codes::AA{Int,3} 
    keys::AA{Tuple{Symbol, Symbol},2}
end

function SynthWorkspaceEclipse(disk::DiskParamsEclipse, lines_number::Int, time_number::Int; ndepths::Integer=100, verbose::Bool=true)
    # allocate the needed memory for synthesis
    lwavgrid = zeros(ndepths)
    rwavgrid = zeros(ndepths)
    allwavs  = zeros(2 * ndepths)
    allints  = zeros(2 * ndepths)
    bist     = zeros(ndepths)
    intt     = zeros(ndepths)
    widt     = zeros(ndepths)

    # allocate the memory for keys, velocities, ld, etc.
    z_cbs = zeros(size(disk.θc)...)
    ϕc = zeros(size(disk.θc)...)
    θc = zeros(size(disk.θc)...)
    μs = zeros(size(disk.θc)..., time_number)
    ld = zeros(size(disk.θc)..., lines_number)
    ext = zeros(size(disk.θc)..., lines_number)
    dA = zeros(size(disk.θc)..., time_number)
    xyz = zeros(size(disk.θc)..., 3) 
    z_rot = zeros(size(disk.θc)...,  time_number)
    ax_codes = zeros(Int, size(disk.θc)..., time_number)
    cbs = zeros(size(disk.θc)...)
    keys = repeat([(:off,:off)], size(disk.θc)...)

    mean_weight_v_no_cb = zeros(length(disk.ϕc), maximum(disk.Nθ), time_number)
    mean_weight_v_earth_orb = zeros(length(disk.ϕc), maximum(disk.Nθ), time_number)
    SP_sun_pos = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_sun_vel = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_bary = fill(Vector{Float64}(undef, 6), disk.Nsubgrid, disk.Nsubgrid)
    pole_vector_grid = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_bary_pos = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_bary_vel = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    OP_bary = fill(Vector{Float64}(undef, 6), disk.Nsubgrid, disk.Nsubgrid)
    projected_velocities_no_cb = zeros(disk.Nsubgrid, disk.Nsubgrid)
    v_scalar_grid = zeros(disk.Nsubgrid, disk.Nsubgrid)
    mu_grid = zeros(disk.Nsubgrid, disk.Nsubgrid)
    v_earth_orb_proj = zeros(disk.Nsubgrid, disk.Nsubgrid)
    distance = zeros(disk.Nsubgrid, disk.Nsubgrid)

    return SynthWorkspaceEclipse(lwavgrid, rwavgrid, allwavs, allints,
                          bist, intt, widt, SP_sun_pos, SP_sun_vel, v_scalar_grid, pole_vector_grid, SP_bary_pos, SP_bary_vel, mu_grid,
                          projected_velocities_no_cb, SP_bary, OP_bary, v_earth_orb_proj, distance, mean_weight_v_no_cb, mean_weight_v_earth_orb,
                          z_cbs, ϕc, θc, μs, ld, ext, dA, xyz, cbs, z_rot, ax_codes, keys)
end
