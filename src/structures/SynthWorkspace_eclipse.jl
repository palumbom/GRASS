struct SynthWorkspaceEclipse{T<:AF}
    lwavgrid::AA{T,1}
    rwavgrid::AA{T,1}
    allwavs::AA{T,1}
    allints::AA{T,1}

    bist::AA{T,1}
    intt::AA{T,1}
    widt::AA{T,1}

    ϕc::AA{T,2}
    θc::AA{T,2}
    μs::AA{T,2}
    ld::AA{T,2}
    dA::AA{T,2}
    wts::AA{T,2}
    xyz::AA{T,3}

    cbs::AA{T,2}
    z_rot::AA{T,2}
    ax_codes::AA{Int,2}
    keys::AA{Tuple{Symbol, Symbol},2}
end

function SynthWorkspaceEclipse(disk::DiskParams; ndepths::Integer=100, verbose::Bool=true)
    # allocate the needed memory for synthesis
    lwavgrid = zeros(ndepths)
    rwavgrid = zeros(ndepths)
    allwavs  = zeros(2 * ndepths)
    allints  = zeros(2 * ndepths)
    bist     = zeros(ndepths)
    intt     = zeros(ndepths)
    widt     = zeros(ndepths)

    # allocate the memory for keys, velocities, ld, etc.
    ϕc = zeros(size(disk.θc)...)
    θc = zeros(size(disk.θc)...)
    μs = zeros(size(disk.θc)...)
    ld = zeros(size(disk.θc)...)
    dA = zeros(size(disk.θc)...)
    xyz = zeros(size(disk.θc)..., 3)
    wts = zeros(size(disk.θc)...)
    z_rot = zeros(size(disk.θc)...)
    ax_codes = zeros(Int, size(disk.θc))

    # allocate additional memory
    cbs = zeros(size(disk.θc)...)
    keys = repeat([(:off,:off)], size(disk.θc)...)

    return SynthWorkspace(lwavgrid, rwavgrid, allwavs, allints,
                          bist, intt, widt, ϕc, θc, μs, ld, dA,
                          wts, xyz, cbs, z_rot, ax_codes, keys)
end
