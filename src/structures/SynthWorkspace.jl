struct SynthWorkspace{T<:AF, N}
    lwavgrid::AA{T,N}
    rwavgrid::AA{T,N}
    allwavs::AA{T,N}
    allints::AA{T,N}
    wavt::AA{T,N}
    bist::AA{T,N}
    dept::AA{T,N}
    widt::AA{T,N}
end

function SynthWorkspace(spec::SpecParams{T}; ndepths::Integer=100) where T
    # get number of dimensions
    if use_gpu
        dims = [length(spec.lines)]
    else
        dims = []
    end

    # allocate the needed memory
    lwavgrid = zeros(ndepths)
    rwavgrid = zeros(ndepths)
    allwavs  = zeros(2 * ndepths)
    allints  = zeros(2 * ndepths)
    wavt     = zeros(ndepths)
    bist     = zeros(ndepths)
    dept     = zeros(ndepths)
    widt     = zeros(ndepths)
    return SynthWorkspace(lwavgrid, rwavgrid, allwavs, allints, wavt, bist, dept, widt)
end
