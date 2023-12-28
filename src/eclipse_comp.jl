function eclipse_compute_quantities!(disk::DiskParams{T}, epoch, obs_long, obs_lat, alt, ϕc::AA{T,2}, θc::AA{T,2},
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

    #query JPL horizons for E, S, M position (km) and velocities (km/s)
    earth_pv = spkssb(399,epoch,"J2000")[1:3] 
    sun_pv = spkssb(10,epoch,"J2000")[1:3] 
    moon_pv = spkssb(301,epoch,"J2000")[1:3] 

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

            # get cartesian coord for each subgrid 
            xyz_sub .= map(x -> sphere_to_cart.(disk.ρs, x...), subgrid) 
            #transform xyz stellar coordinates of grid from sun frame to ICRF
            xyz_sub_bary = map(x -> pxform("IAU_SUN", "J2000", epoch) * x, xyz_sub)
    
            #determine xyz earth coordinates for lat/long of observatory
            EO_earth = pgrrec("EARTH", deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, (earth_radius - earth_radius_pole) / earth_radius)
            #transform xyz earth coordinates of observatory from earth frame to ICRF
            EO_bary = pxform("IAU_EARTH", "J2000", epoch)*EO_earth

            #get vector from barycenter to observatory on Earth's surface
            BO_bary = earth_pv .+ EO_bary
            #get vector from observer to Sun's center 
            OS_bary = BO_bary - sun_pv
            #get vector from observatory on earth's surface to moon center
            OM_bary = moon_pv .- BO_bary
            #get vector from barycenter to each patch on Sun's surface
            BP_bary = map(x -> sun_pv + x, xyz_sub_bary)
            #vectors from observatory on Earth's surface to each patch on Sun's surface
            OP_bary = map(x -> x .- BO_bary, BP_bary)
   

            # calculate mu at each point
            μs_sub .= map(x -> calc_mu(x, OS_bary), xyz_sub_bary) 
            # move to next iteration if patch element is not visible
            all(μs_sub .< zero(T)) && continue


            #determine patches that are blocked by moon 
            #calculate the distance between tile corner and moon
            distance = map(x -> calc_proj_dist(x, OM_bary), OP_bary)

            #get indices for visible patches
            idx1 = μs_sub .> 0.0
            idx3 = (idx1) .& (distance .> atan((moon_radius)/norm(OM_bary))) 

            # assign the mean mu as the mean of visible mus
            μs[i,j] = mean(view(μs_sub, idx3))

            # find xz at mean value of mu and get axis code (i.e., N, E, S, W)
            xyz[i,j,1] = mean(view(getindex.(xyz_sub_bary,1), idx3))
            xyz[i,j,2] = mean(view(getindex.(xyz_sub_bary,2), idx3))
            xyz[i,j,3] = mean(view(getindex.(xyz_sub_bary,3), idx3))
            #ax_codes[i,j] = find_nearest_ax_code(xyz[i,j,1], xyz[i,j,2])

            # calc limb darkening
            # add extinction later on
            #zenith_angle_matrix = rad2deg.(map(x -> calc_proj_dist(x, EO_bary), OP_bary))
            ld_sub .= map(x -> quad_limb_darkening(x), μs_sub)

            # calculate area element of tile
            dA_sub .= map(x -> calc_dA(disk.ρs, getindex(x,1), step(ϕe_sub), step(θe_sub)), subgrid)
            dp_sub .= map((x,y) -> abs(dot(x,y)), OP_bary, xyz_sub_bary) ./ (norm.(OP_bary) .* norm.(xyz_sub_bary))

            # get total projected, visible area of larger tile
            dA_total_proj = dA_sub .* dp_sub

            # copy to workspace
            ld[i,j] = mean(view(ld_sub, idx3))
            dA[i,j] = sum(view(dA_total_proj, idx1))


            # get rotational velocity for location on disk
            z_rot_sub .= map((x,y) -> patch_velocity_los(x..., disk, epoch, y), subgrid, OP_bary)

            wts[i,j] = mean(view(ld_sub .* dA_total_proj, idx3))
            z_rot[i,j] = sum(view(z_rot_sub .* ld_sub, idx3)) ./ sum(view(ld_sub, idx3))
        end
    end
    return nothing
end