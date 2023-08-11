struct GPUAllocs{T1<:AF}
    λs::CuArray{T1,1}

    μs::CuArray{T1,2}
    wts::CuArray{T1,2}
    z_rot::CuArray{T1,2}
    z_cbs::CuArray{T1,2}
    ax_codes::CuArray{Int32,2}

    tloop::CuArray{Int32,2}
    tloop_init::CuArray{Int32,2}
    dat_idx::CuArray{Int32,2}

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
        λs_gpu = CuArray{precision}(spec.lambdas)
        Nθ_gpu = CuArray{Int}(disk.Nθ)
        R_x_gpu = CuArray{precision}(disk.R_x)
        O⃗_gpu = CuArray{precision}(disk.O⃗)
    end

    # allocate memory for synthesis on the GPU
    @cusync begin
        # redshifts, weights, μs
        μs = CUDA.zeros(precision, size(disk.θc))
        wts = CUDA.zeros(precision, size(disk.θc))
        z_rot = CUDA.zeros(precision, size(disk.θc))
        z_cbs = CUDA.zeros(precision, size(disk.θc))
        ax_code = CUDA.zeros(Int32, size(disk.θc))

        # indices
        tloop_gpu = CUDA.zeros(Int32, size(disk.θc))
        tloop_init = CUDA.zeros(Int32, size(disk.θc))
        dat_idx = CUDA.zeros(Int32, size(disk.θc))

        # pre-allocated memory for interpolations
        starmap = CUDA.ones(precision, size(disk.θc)..., Nλ)
        allwavs = CUDA.zeros(precision, size(disk.θc)..., 200)
        allints = CUDA.zeros(precision, size(disk.θc)..., 200)
    end

    return GPUAllocs(λs_gpu, μs, wts, z_rot, z_cbs, ax_code, dat_idx,
                     tloop_gpu, tloop_init, starmap, allwavs, allints)
end
