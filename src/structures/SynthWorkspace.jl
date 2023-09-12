struct SynthWorkspace{T<:AF}
    lwavgrid::AA{T,1}
    rwavgrid::AA{T,1}
    allwavs::AA{T,1}
    allints::AA{T,1}
    bist::AA{T,1}
    intt::AA{T,1}
    widt::AA{T,1}
    μs::AA{T,1}
    cbs::AA{T,1}
    ld::AA{T,1}
    dA::AA{T,1}
    wts::AA{T,1}
    xys::AA{T,2}
    z_rot::AA{T,1}
    ax_codes::AA{Int,1}
    keys::AA{Tuple{Symbol, Symbol},1}
end

function SynthWorkspace(disk::DiskParams; ndepths::Integer=100, verbose::Bool=true)
    # allocate the needed memory for synthesis
    lwavgrid = zeros(ndepths)
    rwavgrid = zeros(ndepths)
    allwavs  = zeros(2 * ndepths)
    allints  = zeros(2 * ndepths)
    bist     = zeros(ndepths)
    intt     = zeros(ndepths)
    widt     = zeros(ndepths)

    # allocate the memory for keys, velocities, ld, etc.
    μs = zeros(size(disk.θc)...)
    ld = zeros(size(disk.θc)...)
    dA = zeros(size(disk.θc)...)
    xyz = zeros(size(disk.θc)..., 3)
    wts = zeros(size(disk.θc)...)
    z_rot = zeros(size(disk.θc)...)
    ax_codes = zeros(Int, size(disk.θc))

    # pre-compute quantities to be re-used
    if verbose
        println("\t>>> Precomputing geometric quantities...")
    end
    precompute_quantities!(disk, μs, ld, dA, xyz, wts, z_rot, ax_codes)

    # get indices with nonzero wts
    idx = μs .> 0.0
    num_nonzero = sum(idx)

    # get arrays of nonzero wts
    μs = μs[idx]
    ld = ld[idx]
    dA = dA[idx]
    wts = wts[idx]
    xyz = xyz[idx, :]
    z_rot = z_rot[idx]
    ax_codes = ax_codes[idx]

    # allocate additional memory
    cbs = zeros(num_nonzero)
    keys = repeat([(:off,:off)], num_nonzero)

    return SynthWorkspace(lwavgrid, rwavgrid, allwavs, allints,
                          bist, intt, widt, μs, cbs, ld, dA,
                          wts, xyz, z_rot, ax_codes, keys)
end
