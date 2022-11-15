struct GPUAllocs{T1<:AF}
    grid::CuArray{T1,1}
    lambdas::CuArray{T1,1}
    tloop::CuArray{Int32,2}
    data_inds::CuArray{Int32,2}
    norm_terms::CuArray{T1,2}
    z_rot::CuArray{T1,2}
    z_cbs::CuArray{T1,2}
    starmap::CuArray{T1,3}
    allwavs::CuArray{T1,3}
    allints::CuArray{T1,3}
end

function GPUAllocs(spec::SpecParams, disk::DiskParams, grid::StepRangeLen; precision::DataType=Float64)
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # move spatial and lambda grid to GPU
    @cusync begin
        grid_gpu = CuArray{precision}(grid)
        lambdas_gpu = CuArray{precision}(spec.lambdas)
    end

    # allocate memory for synthesis on the GPU
    @cusync begin
        # indices, redshifts, and limb darkening
        tloop_gpu = CUDA.zeros(Int32, N, N)
        data_inds = CUDA.zeros(Int32, N, N)
        norm_terms = CUDA.zeros(precision, N, N)
        z_rot = CUDA.zeros(precision, N, N)
        z_cbs = CUDA.zeros(precision, N, N)

        # pre-allocated memory for interpolations
        starmap = CUDA.ones(precision, N, N, Nλ)
        allwavs = CUDA.zeros(precision, N, N, 200)
        allints = CUDA.zeros(precision, N, N, 200)
    end
    return GPUAllocs(grid_gpu, lambdas_gpu, tloop_gpu, data_inds, norm_terms,
                     z_rot, z_cbs, starmap, allwavs, allints)
end
