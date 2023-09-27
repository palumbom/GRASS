function eclipse_compute_quantities!(xyz_planet::AA{T,1}, planet::Planet{T}, # xyz moon and radius instead
    disk::DiskParams{T}, wsp::SynthWorkspace{T}, #observor and memory allocation, wsp - call regular compute quantities inside like rossiter 
    ros_allocs::RossiterAllocs{T}) where T<:AF #memory allocation
# parse out composite type fields
Nsubgrid = disk.Nsubgrid

# allocate memory that wont be needed outside this function
d2_sub = ros_allocs.d2_sub
μs_sub = ros_allocs.μs_sub
ld_sub = ros_allocs.ld_sub
dA_sub = ros_allocs.dA_sub
dp_sub = ros_allocs.dp_sub
xyz_sub = repeat([zeros(3)], Nsubgrid, Nsubgrid)
z_rot_sub = ros_allocs.z_rot_sub
idx1 = ros_allocs.idx1
idx2 = ros_allocs.idx2
idx3 = ros_allocs.idx3

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
θ_l = wsp.θc[i] - dθ/2.0
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
μs_sub .= map(x -> calc_mu(x, disk.O⃗), xyz_sub) #disk.O⃗ varies for each time step - needs to be updated 
if all(μs_sub .<= zero(T))
continue
end

# calculate the distance between subtile center and planet
d2_sub .= map(x -> calc_proj_dist2(x, xyz_planet), xyz_sub)

# if entire course tile visible, use old weights and move on
if all(d2_sub .> planet.radius^2.0) | (xyz_planet[3] < 0.0)
continue
end

# calc limb darkening
ld_sub .= map(x -> quad_limb_darkening(x, disk.u1, disk.u2), μs_sub)

# get rotational velocity for location on disk
z_rot_sub .= map(x -> patch_velocity_los(x..., disk), subgrid)

#2 (this file) 2b: get "propoer" motion between sun and earth - relative velocity is difference between bary velocity of both objects 

# calculate area element of tile
dA_sub .= map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕe_sub), step(θe_sub)), subgrid)
dp_sub .= map(x -> abs(dot(x .- disk.O⃗, x)), xyz_sub) / norm(disk.O⃗)

# get indices for visible patches
idx1 .= μs_sub .> 0.0
idx2 .= d2_sub .> planet.radius^2.0
idx3 .= idx1 .& idx2

# if no patches are visible, set wts, etc. to zero and move on
if all(iszero(idx3))
ros_allocs.μs[i] = 0.0
ros_allocs.ld[i] = 0.0
ros_allocs.dA[i] = 0.0
ros_allocs.wts[i] = 0.0
ros_allocs.z_rot[i] = 0.0
continue
end

# get total projected, visible area of larger tile
dA_total = sum(view(dA_sub, idx3))
dA_total_proj = sum(view(dA_sub .* dp_sub, idx3))

# set limb darkening as mean of visible patches
ros_allocs.μs[i] = mean(view(μs_sub, idx1))
ros_allocs.ld[i] = mean(view(ld_sub, idx3))
ros_allocs.dA[i] = dA_total_proj

ros_allocs.wts[i] = mean(view(ld_sub .* dA_total_proj, idx3))
ros_allocs.z_rot[i] = sum(view(z_rot_sub .* ld_sub, idx3)) ./ sum(view(ld_sub, idx3))
end
return nothing
end
