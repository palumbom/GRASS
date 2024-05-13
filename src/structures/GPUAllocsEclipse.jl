struct GPUAllocsEclipse{T1<:AF}
    λs::CuArray{T1,1}
    prof::CuArray{T1,1}
    flux::CuArray{T1,2}

    μs::CuArray{T1,2}
    dA::CuArray{T1,2}
    ld::CuArray{T1,3}
    ext::CuArray{T1,3}
    z_rot::CuArray{T1,2}
    z_cbs::CuArray{T1,2}

    tloop::CuArray{Int32,2}
    dat_idx::CuArray{Int32,2}
    ax_codes::CuArray{Int32,2}

    allwavs::CuArray{T1,3}
    allints::CuArray{T1,3}
end

function GPUAllocsEclipse(spec::SpecParams, disk::DiskParamsEclipse, lines_number::Int; precision::DataType=Float64, verbose::Bool=true)
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    @cusync begin
        λs_gpu = CuArray{precision}(spec.lambdas)
        prof_gpu = CUDA.zeros(precision, Nλ)
        flux_gpu = CUDA.ones(precision, Nλ, Nt)
    end

    # allocate memory for precomputation
    Nϕ = disk.N
    Nθ_max = maximum(disk.Nθ)

    # allocate memory for pre-computations
    @cusync begin
        μs = CUDA.zeros(precision, Nϕ, Nθ_max)
        dA = CUDA.zeros(precision, Nϕ, Nθ_max)
        ld = CUDA.zeros(precision, Nϕ, Nθ_max, lines_number)
        ext = CUDA.zeros(precision, Nϕ, Nθ_max, lines_number)
        z_rot = CUDA.zeros(precision, Nϕ, Nθ_max)
        ax_code = CUDA.zeros(Int32, Nϕ, Nθ_max)
    end

    # allocate memory for convective blueshifts
    @cusync z_cbs = CUDA.zeros(precision, Nϕ, Nθ_max)

    # allocate memory for indices
    @cusync begin
        tloop_gpu = CUDA.zeros(Int32, Nϕ, Nθ_max)
        dat_idx = CUDA.zeros(Int32, Nϕ, Nθ_max)
    end

    # allocated memory for synthesis
    @cusync begin
        allwavs = CUDA.zeros(precision, Nϕ, Nθ_max, 200)
        allints = CUDA.zeros(precision, Nϕ, Nθ_max, 200)
    end

    return GPUAllocsEclipse(λs_gpu, prof_gpu, flux_gpu, μs, dA, ld, ext, z_rot, z_cbs, tloop_gpu,
                     dat_idx, ax_code, allwavs, allints)
end