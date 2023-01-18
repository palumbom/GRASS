struct SynthWorkspace{T<:AF, N}
    lwavgrid::AA{T,N}
    rwavgrid::AA{T,N}
    allwavs::AA{T,N}
    allints::AA{T,N}
    bist::AA{T,N}
    intt::AA{T,N}
    widt::AA{T,N}
end

function SynthWorkspace(;ndepths::Integer=100) where T<:AF
    # allocate the needed memory
    lwavgrid = zeros(T, ndepths)
    rwavgrid = zeros(T, ndepths)
    allwavs  = zeros(T, 2 * ndepths)
    allints  = zeros(T, 2 * ndepths)
    bist     = zeros(T, ndepths)
    intt     = zeros(T, ndepths)
    widt     = zeros(T, ndepths)
    return SynthWorkspace(lwavgrid, rwavgrid, allwavs, allints, bist, intt, widt)
end
