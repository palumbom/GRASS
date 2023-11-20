# function eclipse_compute_quantities!(xyz_planet::AA{T,1}, r_moon, # xyz moon and radius instead
#                                      disk::DiskParams{T}, wsp::SynthWorkspace{T}, #observor and memory allocation, wsp - call regular compute quantities inside like rossiter
#                                      ros_allocs::RossiterAllocs{T}) where T<:AF #memory allocation
#     # parse out composite type fields
#     Nsubgrid = disk.Nsubgrid

#     # allocate memory that wont be needed outside this function
#     d2_sub = ros_allocs.d2_sub
#     μs_sub = ros_allocs.μs_sub
#     ld_sub = ros_allocs.ld_sub
#     dA_sub = ros_allocs.dA_sub
#     dp_sub = ros_allocs.dp_sub
#     xyz_sub = repeat([zeros(3)], Nsubgrid, Nsubgrid)
#     z_rot_sub = ros_allocs.z_rot_sub
#     idx1 = ros_allocs.idx1
#     idx2 = ros_allocs.idx2
#     idx3 = ros_allocs.idx3

#     # loop over disk positions
#     for i in eachindex(wsp.ϕc)
#         # get the number of theta tiles needed for the latitude tiles
#         Nθ = get_Nθ(wsp.ϕc[i], step(disk.ϕe))

#         # get the latitude and longitude increments
#         dϕ = step(disk.ϕe)
#         dθ = deg2rad(360.0) / Nθ

#         # get edges of large tile
#         ϕ_l = wsp.ϕc[i] - dϕ/2.0
#         ϕ_r = wsp.ϕc[i] + dϕ/2.0
#         θ_l = wsp.θc[i] - dθ/2.0
#         θ_r = wsp.θc[i] + dθ/2.0

#         # subdivide the tile
#         ϕe_sub = range(ϕ_l, ϕ_r, length=Nsubgrid+1)
#         θe_sub = range(θ_l, θ_r, length=Nsubgrid+1)
#         ϕc_sub = get_grid_centers(ϕe_sub)
#         θc_sub = get_grid_centers(θe_sub)
#         subgrid = Iterators.product(ϕc_sub, θc_sub)

#         # get cartesian coord for each subgrid and rotate by rot. matrix
#         xyz_sub .= map(x -> sphere_to_cart.(disk.ρs, x...), subgrid)
#         xyz_sub .= map(x -> disk.R_x * x, xyz_sub)

#         # calculate mu at each point
#         μs_sub .= map(x -> calc_mu(x, disk.O⃗), xyz_sub) #disk.O⃗ varies for each time step - needs to be updated
#         if all(μs_sub .<= zero(T))
#             continue
#         end

#         # calculate the distance between subtile center and planet
#         d2_sub .= map(x -> calc_proj_dist2(x, xyz_planet), xyz_sub)

#         # if entire course tile visible, use old weights and move on
#         if all(d2_sub .> planet.radius^2.0) | (xyz_planet[3] < 0.0)
#             continue
#         end

#         # calc limb darkening
#         ld_sub .= map(x -> quad_limb_darkening(x, disk.u1, disk.u2), μs_sub)

#         # get rotational velocity for location on disk
#         z_rot_sub .= map(x -> patch_velocity_los(x..., disk), subgrid)

#         #2 (this file) 2b: get "propoer" motion between sun and earth - relative velocity is difference between bary velocity of both objects

#         # calculate area element of tile
#         dA_sub .= map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕe_sub), step(θe_sub)), subgrid)
#         dp_sub .= map(x -> abs(dot(x .- disk.O⃗, x)), xyz_sub) / norm(disk.O⃗)

#         # get indices for visible patches
#         idx1 .= μs_sub .> 0.0
#         idx2 .= d2_sub .> planet.radius^2.0
#         idx3 .= idx1 .& idx2

#         # if no patches are visible, set wts, etc. to zero and move on
#         if all(iszero(idx3))
#             ros_allocs.μs[i] = 0.0
#             ros_allocs.ld[i] = 0.0
#             ros_allocs.dA[i] = 0.0
#             ros_allocs.wts[i] = 0.0
#             ros_allocs.z_rot[i] = 0.0
#             continue
#         end

#         # get total projected, visible area of larger tile
#         dA_total = sum(view(dA_sub, idx3))
#         dA_total_proj = sum(view(dA_sub .* dp_sub, idx3))

#         # set limb darkening as mean of visible patches
#         ros_allocs.μs[i] = mean(view(μs_sub, idx1))
#         ros_allocs.ld[i] = mean(view(ld_sub, idx3))
#         ros_allocs.dA[i] = dA_total_proj

#         ros_allocs.wts[i] = mean(view(ld_sub .* dA_total_proj, idx3))
#         ros_allocs.z_rot[i] = sum(view(z_rot_sub .* ld_sub, idx3)) ./ sum(view(ld_sub, idx3))
#     end
#     return nothing
# end

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
