struct GPUAllocs{T1<:AF}
    ϕc::CuArray{T1,1}
    θc::CuArray{T1,1}
    R_θ::CuArray{T1,2}
    O⃗::CuArray{T1,1}
    λs::CuArray{T1,1}
    μs::CuArray{T1,2}
    vec1::CuArray{T1,3}
    vec2::CuArray{T1,3}
    vec3::CuArray{T1,3}
    tloop::CuArray{Int32,2}
    dat_idx::CuArray{Int32,2}
    weights::CuArray{T1,2}
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
        # vector storage space
        vec1 = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc), 3)
        vec2 = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc), 3)
        vec3 = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc), 3)

        # indices
        tloop_gpu = CUDA.zeros(Int32, length(disk.ϕc), length(disk.θc))
        dat_idx = CUDA.zeros(Int32, length(disk.ϕc), length(disk.θc))

        # redshifts, weights, μs
        μs_gpu = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))
        weights = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))
        z_rot = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))
        z_cbs = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc))

        # pre-allocated memory for interpolations
        starmap = CUDA.ones(precision, length(disk.ϕc), length(disk.θc), Nλ)
        allwavs = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc), 200)
        allints = CUDA.zeros(precision, length(disk.ϕc), length(disk.θc), 200)
    end
    return GPUAllocs(ϕ_gpu, θ_gpu, R_θ_gpu, O⃗_gpu, λs_gpu, μs_gpu,
                     vec1, vec2, vec3, tloop_gpu, dat_idx, weights,
                     z_rot, z_cbs, starmap, allwavs, allints)
end
