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
    lwavgrid = ArrayType(zeros(ndepths), dims...)
    rwavgrid = ArrayType(zeros(ndepths), dims...)
    allwavs  = ArrayType(zeros(2 * ndepths), dims...)
    allints  = ArrayType(zeros(2 * ndepths), dims...)
    wavt     = ArrayType(zeros(ndepths), dims...)
    bist     = ArrayType(zeros(ndepths), dims...)
    dept     = ArrayType(zeros(ndepths), dims...)
    widt     = ArrayType(zeros(ndepths), dims...)
    return SynthWorkspace(lwavgrid, rwavgrid, allwavs, allints, wavt, bist, dept, widt)
end
