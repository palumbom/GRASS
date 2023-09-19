function eclipse_compute_quantities!(disk::DiskParams{T}, ϕc::AA{T,2}, θc::AA{T,2},
    μs::AA{T,2}, ld::AA{T,2}, dA::AA{T,2},
    xyz::AA{T,3}, wts::AA{T,2}, z_rot::AA{T,2},
    ax_codes::AA{Int64, 2}) where T<:AF
# parse out composite type fields
Nsubgrid = disk.Nsubgrid

# allocate memory that wont be needed outside this function
μs_sub = zeros(Nsubgrid, Nsubgrid)
ld_sub = zeros(Nsubgrid, Nsubgrid)
dA_sub = zeros(Nsubgrid, Nsubgrid)
dp_sub = zeros(Nsubgrid, Nsubgrid)
xyz_sub = repeat([zeros(3)], Nsubgrid, Nsubgrid)
z_rot_sub = zeros(Nsubgrid, Nsubgrid)
idx = BitMatrix(undef, size(μs_sub))

# loop over disk positions
for i in eachindex(disk.ϕc)
for j in 1:disk.Nθ[i]
# save the tile position
ϕc[i,j] = disk.ϕc[i]
θc[i,j] = disk.θc[i,j]

# subdivide the tile
ϕe_sub = range(disk.ϕe[i], disk.ϕe[i+1], length=Nsubgrid+1)
θe_sub = range(disk.θe[i,j], disk.θe[i,j+1], length=Nsubgrid+1)
ϕc_sub = get_grid_centers(ϕe_sub)
θc_sub = get_grid_centers(θe_sub)
subgrid = Iterators.product(ϕc_sub, θc_sub)

# get cartesian coord for each subgrid and rotate by rot. matrix
xyz_sub .= map(x -> sphere_to_cart.(disk.ρs, x...), subgrid)
xyz_sub .= map(x -> disk.R_x * x, xyz_sub)

# calculate mu at each point
μs_sub .= map(x -> calc_mu(x, disk.O⃗), xyz_sub)

# move to next iteration if patch element is not visible
all(μs_sub .<= zero(T)) && continue

# assign the mean mu as the mean of visible mus
idx .= μs_sub .> 0.0
μs[i,j] = mean(view(μs_sub, idx))

# find xz at mean value of mu and get axis code (i.e., N, E, S, W)
xyz[i,j,1] = mean(view(getindex.(xyz_sub,1), idx))
xyz[i,j,2] = mean(view(getindex.(xyz_sub,2), idx))
xyz[i,j,3] = mean(view(getindex.(xyz_sub,3), idx))
ax_codes[i,j] = find_nearest_ax_code(xyz[i,j,1], xyz[i,j,2])

# calc limb darkening
ld_sub .= map(x -> quad_limb_darkening(x, disk.u1, disk.u2), μs_sub)

# get rotational velocity for location on disk
z_rot_sub .= map(x -> patch_velocity_los(x..., disk), subgrid)

# calculate area element of tile
dA_sub .= map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕe_sub), step(θe_sub)), subgrid)
dp_sub .= map(x -> abs(dot(x .- disk.O⃗, x)), xyz_sub) ./ norm(disk.O⃗)

# get total projected, visible area of larger tile
dA_total = sum(view(dA_sub, idx))
dA_total_proj = sum(view(dA_sub .* dp_sub, idx))

# copy to workspace
ld[i,j] = mean(view(ld_sub, idx))
dA[i,j] = dA_total_proj

wts[i,j] = mean(view(ld_sub .* dA_total_proj, idx))
z_rot[i,j] = sum(view(z_rot_sub .* ld_sub, idx)) ./ sum(view(ld_sub, idx))
end
end
return nothing
end