struct RossiterAllocs{T<:AF}
    ϕc::AA{T,1}
    θc::AA{T,1}
    μs::AA{T,1}
    ld::AA{T,1}
    dA::AA{T,1}
    wts::AA{T,1}
    xyz::AA{T,2}
end

function RossiterAllocs(wsp::SynthWorkspace{T}) where T<:AF
    return RossiterAllocs(deepcopy(wsp.ϕc),
                          deepcopy(wsp.θc),
                          deepcopy(wsp.μs),
                          deepcopy(wsp.ld),
                          deepcopy(wsp.dA),
                          deepcopy(wsp.wts),
                          deepcopy(wsp.xyz))
end
