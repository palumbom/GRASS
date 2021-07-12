struct SynthWorkspace{T<:AF}
    lwavgrid::AA{T,1}
    rwavgrid::AA{T,1}
    allwavs::AA{T,1}
    allints::AA{T,1}
    wavt::AA{T,1}
    bist::AA{T,1}
    dept::AA{T,1}
    widt::AA{T,1}
end

function SynthWorkspace(; ndepths::Integer=100)
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
