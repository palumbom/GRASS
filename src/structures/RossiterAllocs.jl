 struct RossiterAllocs{T<:AF}
    μs::AA{T,1}
    ld::AA{T,1}
    dA::AA{T,1}
    wts::AA{T,1}
    z_rot::AA{T,1}

    d2_sub::AA{T,2}
    μs_sub::AA{T,2}
    ld_sub::AA{T,2}
    dA_sub::AA{T,2}
    dp_sub::AA{T,2}
    # xyz_sub::AA{T,2}
    z_rot_sub::AA{T,2}

    idx1::BitArray{2}
    idx2::BitArray{2}
    idx3::BitArray{2}
end

function RossiterAllocs(wsp::SynthWorkspace{T}, disk::DiskParams{T}) where T<:AF
    # get number of grid cells needed
    Nsubgrid = disk.Nsubgrid

    # allocate memory for computations on subgrid
    d2_sub = zeros(Nsubgrid, Nsubgrid)
    μs_sub = zeros(Nsubgrid, Nsubgrid)
    ld_sub = zeros(Nsubgrid, Nsubgrid)
    dA_sub = zeros(Nsubgrid, Nsubgrid)
    dp_sub = zeros(Nsubgrid, Nsubgrid)
    xyz_sub = repeat([zeros(3)], Nsubgrid, Nsubgrid)
    z_rot_sub = zeros(Nsubgrid, Nsubgrid)
    idx1 = BitMatrix(undef, size(μs_sub))
    idx2 = BitMatrix(undef, size(μs_sub))
    idx3 = BitMatrix(undef, size(μs_sub))

    # copy things from workspace
    μs = copy(wsp.μs)
    ld = copy(wsp.ld)
    dA = copy(wsp.dA)
    wts = copy(wsp.wts)
    z_rot = copy(wsp.z_rot)

    return RossiterAllocs(μs, ld, dA, wts, z_rot,
                          d2_sub, μs_sub, ld_sub,
                          dA_sub, dp_sub, z_rot_sub,
                          idx1, idx2, idx3)
end

struct RossiterAllocsGPU{T<:AF}
    μs::CuArray{T,1}
    wts::CuArray{T,1}
    z_rot::CuArray{T,1}
end

function RossiterAllocsGPU(gpu_allocs::GPUAllocs{T}) where T<:AF
    @cusync begin
        μs = copy(gpu_allocs.μs)
        wts = copy(gpu_allocs.wts)
        z_rot = copy(gpu_allocs.z_rot)
    end
    return RossiterAllocsGPU(μs, wts, z_rot)
end
