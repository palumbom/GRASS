struct GeoWorkspaceEclipse{T<:AF}
    dA_total_proj_mean::AA{T,3}
    mean_intensity::AA{T,3}
    mean_weight_v_no_cb::AA{T,3}
    mean_weight_v_earth_orb::AA{T,3}
    zenith_mean::AA{T,3}

    pole_vector_grid::Matrix{Vector{Float64}}
    SP_sun_pos::Matrix{Vector{Float64}}
    SP_sun_vel::Matrix{Vector{Float64}}
    SP_bary::Matrix{Vector{Float64}}
    SP_bary_pos::Matrix{Vector{Float64}}
    SP_bary_vel::Matrix{Vector{Float64}}
    OP_bary::Matrix{Vector{Float64}}

    mu_grid::AA{T,2}
    projected_velocities_no_cb::AA{T,2}
    distance::AA{T,2}
    v_scalar_grid::AA{T,2}
    v_earth_orb_proj::Matrix{Float64}
end

function GeoWorkspaceEclipse(disk::DiskParamsEclipse, lines_number::Int, time_number::Int)
    # allocate memory that wont be needed outside this function
    dA_total_proj_mean = zeros(length(disk.ϕc), maximum(disk.Nθ), time_number)
    mean_intensity = zeros(length(disk.ϕc), maximum(disk.Nθ), lines_number)
    mean_weight_v_no_cb = zeros(length(disk.ϕc), maximum(disk.Nθ), time_number)
    mean_weight_v_earth_orb = zeros(length(disk.ϕc), maximum(disk.Nθ), time_number)
    zenith_mean = zeros(length(disk.ϕc), maximum(disk.Nθ), time_number)
    
    # vectors
    pole_vector_grid = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_sun_pos = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_sun_vel = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_bary = fill(Vector{Float64}(undef, 6), disk.Nsubgrid, disk.Nsubgrid)
    SP_bary_pos = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    SP_bary_vel = fill(Vector{Float64}(undef, 3), disk.Nsubgrid, disk.Nsubgrid)
    OP_bary = fill(Vector{Float64}(undef, 6), disk.Nsubgrid, disk.Nsubgrid)
    # scalars
    mu_grid = zeros(disk.Nsubgrid, disk.Nsubgrid)
    projected_velocities_no_cb = zeros(disk.Nsubgrid, disk.Nsubgrid)
    distance = zeros(disk.Nsubgrid, disk.Nsubgrid)
    v_scalar_grid = zeros(disk.Nsubgrid, disk.Nsubgrid)
    v_earth_orb_proj = zeros(disk.Nsubgrid, disk.Nsubgrid)

    return GeoWorkspaceEclipse(dA_total_proj_mean, mean_intensity, mean_weight_v_no_cb, mean_weight_v_earth_orb,
    zenith_mean, pole_vector_grid, SP_sun_pos, SP_sun_vel, SP_bary, SP_bary_pos, SP_bary_vel, OP_bary, mu_grid,
    projected_velocities_no_cb, distance, v_scalar_grid, v_earth_orb_proj)
end