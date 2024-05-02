struct EclipseAllocs{T<:AF}
    lwavgrid::CuArray{T,1}
    rwavgrid::CuArray{T,1}
    allwavs::CuArray{T,1}
    allints::CuArray{T,1}

    bist::CuArray{T,1}
    intt::CuArray{T,1}
    widt::CuArray{T,1}

    ϕc::CuArray{T,2}
    θc::CuArray{T,2}
    μs::CuArray{T,2}
    ld::CuArray{T,3}
    ext::CuArray{T,3}
    dA::CuArray{T,2}
    wts::CuArray{T,3}
    xyz::CuArray{T,3}

    cbs::CuArray{T,2}
    z_rot::CuArray{T,3}
    ax_codes::CuArray{Int,2}
    keys::CuArray{Tuple{Symbol, Symbol},2}

    dA_total_proj_mean::CuArray{T,2}
    mean_intensity::CuArray{T,2}
    mean_weight_v_no_cb::CuArray{T,2}
    mean_weight_v_earth_orb::CuArray{T,2}

    pole_vector_grid::CuArray{CuArray{Float64}}
    SP_sun_pos::CuArray{CuArray{Float64}}
    SP_sun_vel::CuArray{CuArray{Float64}}
    SP_bary::CuArray{CuArray{Float64}}
    SP_bary_pos::CuArray{CuArray{Float64}}
    SP_bary_vel::CuArray{CuArray{Float64}}
    OP_bary::CuArray{CuArray{Float64}}

    mu_grid::CuArray{T,2}
    projected_velocities_no_cb::CuArray{T,2}
    distance::CuArray{T,2}
    v_scalar_grid::CuArray{T,2}
    v_earth_orb_proj::CuArray{Float64}
end

function EclipseAllocs(wsp::SynthWorkspaceEclipse{T1}, mem::GeoWorkspaceEclipse{T2}) where {T1<:AF, T2<:AF}
    lwavgrid = CuArray(wsp.lwavgrid)
    rwavgrid = CuArray(wsp.rwavgrid)
    allwavs = CuArray(wsp.allwavs)
    allints = CuArray(wsp.allints)

    bist = CuArray(wsp.bist)
    intt = CuArray(wsp.intt)
    widt = CuArray(wsp.widt)

    ϕc = CuArray(wsp.ϕc)
    θc = CuArray(wsp.θc)
    μs = CuArray(wsp.μs)
    ld = CuArray(wsp.ld)
    ext = CuArray(wsp.ext)
    dA = CuArray(wsp.dA)
    wts = CuArray(wsp.wts)
    xyz = CuArray(wsp.xyz)

    cbs = CuArray(wsp.cbs)
    z_rot = CuArray(wsp.z_rot)
    ax_codes = CuArray(wsp.ax_codes)
    keys = CuArray(wsp.keys)

    dA_total_proj_mean = CuArray(mem.dA_total_proj_mean)
    mean_intensity = CuArray(mem.mean_intensity)
    mean_weight_v_no_cb = CuArray(mem.mean_weight_v_no_cb)
    mean_weight_v_earth_orb = CuArray(mem.mean_weight_v_earth_orb)

    pole_vector_grid = CuArray(mem.pole_vector_grid)
    SP_sun_pos = CuArray(mem.SP_sun_pos)
    SP_sun_vel = CuArray(mem.SP_sun_vel)
    SP_bary = CuArray(mem.SP_bary)
    SP_bary_pos = CuArray(mem.SP_bary_pos)
    SP_bary_vel = CuArray(mem.SP_bary_vel)
    OP_bary = CuArray(mem.OP_bary)

    mu_grid = CuArray(mem.mu_grid)
    projected_velocities_no_cb = CuArray(mem.projected_velocities_no_cb)
    distance = CuArray(mem.distance)
    v_scalar_grid = CuArray(mem.v_scalar_grid)
    v_earth_orb_proj = CuArray(mem.v_earth_orb_proj)

    return EclipseAllocs(lwavgrid, rwavgrid, allwavs, allints, bist, intt, widt, ϕc, θc, μs, ld, \
                            ext, dA, wts, xyz, cbs, z_rot, ax_codes, keys, dA_total_proj_mean, \
                            mean_intensity, mean_weight_v_no_cb, mean_weight_v_earth_orb, \
                            pole_vector_grid, SP_sun_pos, SP_sun_vel, SP_bary, SP_bary_pos, \
                            SP_bary_vel, OP_bary, mu_grid, projected_velocities_no_cb, distance, \
                            v_scalar_grid, v_earth_orb_proj)
end