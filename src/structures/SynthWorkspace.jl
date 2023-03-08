struct SynthWorkspace{T<:AF, N}
    lwavgrid::AA{T,N}
    rwavgrid::AA{T,N}
    allwavs::AA{T,N}
    allints::AA{T,N}
    bist::AA{T,N}
    intt::AA{T,N}
    widt::AA{T,N}
end

function SynthWorkspace(;ndepths::Integer=100)
    # allocate the needed memory
    lwavgrid = zeros(ndepths)
    rwavgrid = zeros(ndepths)
    allwavs  = zeros(2 * ndepths)
    allints  = zeros(2 * ndepths)
    bist     = zeros(ndepths)
    intt     = zeros(ndepths)
    widt     = zeros(ndepths)
    return SynthWorkspace(lwavgrid, rwavgrid, allwavs, allints, bist, intt, widt)
end
