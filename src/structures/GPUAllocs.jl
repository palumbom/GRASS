struct GPUAllocs{T1<:AF}
    λs::CuArray{T1,1}

    μs::CuArray{T1,1}
    wts::CuArray{T1,1}
    z_rot::CuArray{T1,1}
    z_cbs::CuArray{T1,1}
    ax_codes::CuArray{Int32,1}

    tloop::CuArray{Int32,1}
    tloop_init::CuArray{Int32,1}
    dat_idx::CuArray{Int32,1}

    starmap::CuArray{T1,2}
    allwavs::CuArray{T1,2}
    allints::CuArray{T1,2}
end

function GPUAllocs(spec::SpecParams, disk::DiskParams; precision::DataType=Float64, verbose::Bool=true)
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

    # pre-compute quantities to be re-used
    if verbose
        println("\t>>> Precomputing geometric quantities...")
    end

    # allocate memory for precomputation
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)
    num_tiles = Nϕ * Nθ_max

    # allocate memory for pre-computations
    @cusync begin
        μs = CUDA.zeros(precision, Nϕ, Nθ_max)
        wts = CUDA.zeros(precision, Nϕ, Nθ_max)
        z_rot = CUDA.zeros(precision, Nϕ, Nθ_max)
        ax_code = CUDA.zeros(Int32, Nϕ, Nθ_max)
    end

    # perform the pre-computations
    precompute_quantities_gpu!(disk, μs, wts, z_rot, ax_code)

    # reshape to a vector
    @cusync μs = CUDA.reshape(μs, num_tiles)
    @cusync wts = CUDA.reshape(wts, num_tiles)
    @cusync z_rot = CUDA.reshape(z_rot, num_tiles)
    @cusync ax_code = CUDA.reshape(ax_code, num_tiles)

    # get indices with nonzero wts
    @cusync idx = wts .!= 0.0
    @cusync num_nonzero = CUDA.sum(idx)

    # filter nonzero weights
    μs = filter_array_gpu(μs, idx)
    wts = filter_array_gpu(wts, idx)
    z_rot = filter_array_gpu(z_rot, idx)
    ax_code = filter_array_gpu(ax_code, idx)

    # allocate memory for convective blueshifts
    @cusync z_cbs = CUDA.zeros(precision, num_nonzero)

    # allocate memory for indices
    @cusync begin
        tloop_gpu = CUDA.zeros(Int32, num_nonzero)
        tloop_init = CUDA.zeros(Int32, num_nonzero)
        dat_idx = CUDA.zeros(Int32, num_nonzero)
    end

    # allocated memory for synthesis
    @cusync begin
        starmap = CUDA.ones(precision, num_nonzero, Nλ)
        allwavs = CUDA.zeros(precision, num_nonzero, 200)
        allints = CUDA.zeros(precision, num_nonzero, 200)
    end

    return GPUAllocs(λs_gpu, μs, wts, z_rot, z_cbs, ax_code, dat_idx,
                     tloop_gpu, tloop_init, starmap, allwavs, allints)
end

function filter_array_gpu(input::CuArray{T,1}, pred::CuArray{Bool, 1}) where T
    # initialize counter
    counter = 0

    # allocate memory for output
    @cusync output = CUDA.zeros(T, CUDA.sum(pred))

    # launch the kernel
    @cusync @cuda filter_array_gpu!(output, input, pred, counter)

    # queue CUDA to gc
    @cusync CUDA.unsafe_free!(input)
    GC.gc(true)
    return output
end
