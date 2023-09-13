 struct RossiterAllocs{T<:AF}
    μs::AA{T,1}
    ld::AA{T,1}
    dA::AA{T,1}
    wts::AA{T,1}
    z_rot::AA{T,1}
end

function RossiterAllocs(wsp::SynthWorkspace{T}) where T<:AF
    return RossiterAllocs(deepcopy(wsp.μs),
                          deepcopy(wsp.ld),
                          deepcopy(wsp.dA),
                          deepcopy(wsp.wts),
                          deepcopy(wsp.z_rot))
end

struct RossiterAllocsGPU{T<:AF}
    μs::CuArray{T,1}
    wts::CuArray{T,1}
    z_rot::CuArray{T,1}
end

function RossiterAllocsGPU(gpu_allocs::GPUAllocs{T}) where T<:AF
    return RossiterAllocsGPU(deepcopy(gpu_allocs.μs),
                          deepcopy(gpu_allocs.wts),
                          deepcopy(gpu_allocs.z_rot))
end
