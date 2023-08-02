struct GPUAllocs{T1<:AF}
    ϕc::CuArray{T1,1}
    θc::CuArray{T1,1}
    R_θ::CuArray{T1,2}
    O⃗::CuArray{T1,1}
    λs::CuArray{T1,1}
    μs::CuArray{T1,2}
    tloop::CuArray{Int32,2}
    data_inds::CuArray{Int32,2}
    norm_terms::CuArray{T1,2}
    z_rot::CuArray{T1,2}
    z_cbs::CuArray{T1,2}
    starmap::CuArray{T1,3}
    allwavs::CuArray{T1,3}
    allints::CuArray{T1,3}
end

function GPUAllocs(spec::SpecParams, disk::DiskParams; precision::DataType=Float64)
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # move disk + geometry information to gpu
    @cusync begin
        ϕ_gpu = CuArray{precision}(disk.ϕc)
        θ_gpu = CuArray{precision}(disk.θc)
        R_θ_gpu = CuArray{precision}(disk.R_θ)
        O⃗_gpu = CuArray{precision}(disk.O⃗)
        λs_gpu = CuArray{precision}(spec.lambdas)
    end

    # allocate memory for synthesis on the GPU
    @cusync begin
        # indices, redshifts, and limb darkening
        tloop_gpu = CUDA.zeros(Int32, length(disk.ϕc), length(disk.θc))
        data_inds = CUDA.zeros(Int32, length(disk.ϕc), length(disk.θc))
        μs_gpu = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))
        norm_terms = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))
        z_rot = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))
        z_cbs = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))

        # pre-allocated memory for interpolations
        starmap = CUDA.ones(precision, length(disk.ϕc), length(disk.θc), Nλ)
        allwavs = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc), 200)
        allints = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc), 200)
    end
    return GPUAllocs(ϕ_gpu, θ_gpu, R_θ_gpu, O⃗_gpu, λs_gpu, μs_gpu,
                     tloop_gpu, data_inds, norm_terms, z_rot,
                     z_cbs, starmap, allwavs, allints)
end
