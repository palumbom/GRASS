struct GPUAllocs{T1<:AF}
    λs::CuArray{T1,1}
    prof::CuArray{T1,1}

    ϕc::CuArray{T1,1}
    θc::CuArray{T1,1}
    μs::CuArray{T1,1}
    xx::CuArray{T1,1}
    yy::CuArray{T1,1}
    zz::CuArray{T1,1}

    wts::CuArray{T1,1}
    z_rot::CuArray{T1,1}
    z_cbs::CuArray{T1,1}
    ax_codes::CuArray{Int32,1}

    tloop::CuArray{Int32,1}
    tloop_init::CuArray{Int32,1}
    dat_idx::CuArray{Int32,1}

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
        Nθ_gpu = CuArray{Int}(disk.Nθ)
        R_x_gpu = CuArray{precision}(disk.R_x)
        O⃗_gpu = CuArray{precision}(disk.O⃗)
    end

    @cusync begin
        λs_gpu = CuArray{precision}(spec.lambdas)
        prof_gpu = CUDA.zeros(precision, Nλ)
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
        ϕc = CUDA.zeros(precision, Nϕ, Nθ_max)
        θc = CUDA.zeros(precision, Nϕ, Nθ_max)
        μs = CUDA.zeros(precision, Nϕ, Nθ_max)
        xx = CUDA.zeros(precision, Nϕ, Nθ_max)
        yy = CUDA.zeros(precision, Nϕ, Nθ_max)
        zz = CUDA.zeros(precision, Nϕ, Nθ_max)
        wts = CUDA.zeros(precision, Nϕ, Nθ_max)
        z_rot = CUDA.zeros(precision, Nϕ, Nθ_max)
        ax_code = CUDA.zeros(Int32, Nϕ, Nθ_max)
    end

    # perform the pre-computations
    precompute_quantities_gpu!(disk, ϕc, θc, xx, yy, zz, μs, wts, z_rot, ax_code)

    # reshape to a vector
    @cusync ϕc = CUDA.reshape(ϕc, num_tiles)
    @cusync θc = CUDA.reshape(θc, num_tiles)
    @cusync μs = CUDA.reshape(μs, num_tiles)
    @cusync xx = CUDA.reshape(xx, num_tiles)
    @cusync yy = CUDA.reshape(yy, num_tiles)
    @cusync zz = CUDA.reshape(zz, num_tiles)
    @cusync wts = CUDA.reshape(wts, num_tiles)
    @cusync z_rot = CUDA.reshape(z_rot, num_tiles)
    @cusync ax_code = CUDA.reshape(ax_code, num_tiles)

    # get indices with nonzero wts
    @cusync idx = μs .> 0.0
    @cusync num_nonzero = CUDA.sum(idx)

    @cusync begin
        ϕc_new = CUDA.zeros(precision, CUDA.sum(idx))
        θc_new = CUDA.zeros(precision, CUDA.sum(idx))
        μs_new = CUDA.zeros(precision, CUDA.sum(idx))
        xx_new = CUDA.zeros(precision, CUDA.sum(idx))
        yy_new = CUDA.zeros(precision, CUDA.sum(idx))
        zz_new = CUDA.zeros(precision, CUDA.sum(idx))
        wts_new = CUDA.zeros(precision, CUDA.sum(idx))
        z_rot_new = CUDA.zeros(precision, CUDA.sum(idx))
        ax_code_new = CUDA.zeros(Int32, CUDA.sum(idx))
    end

    # launch the kernel
    @cusync @cuda filter_array_gpu!(ϕc_new, ϕc, idx, 0)
    @cusync @cuda filter_array_gpu!(θc_new, θc, idx, 0)
    @cusync @cuda filter_array_gpu!(μs_new, μs, idx, 0)
    @cusync @cuda filter_array_gpu!(xx_new, xx, idx, 0)
    @cusync @cuda filter_array_gpu!(yy_new, yy, idx, 0)
    @cusync @cuda filter_array_gpu!(zz_new, zz, idx, 0)
    @cusync @cuda filter_array_gpu!(wts_new, wts, idx, 0)
    @cusync @cuda filter_array_gpu!(z_rot_new, z_rot, idx, 0)
    @cusync @cuda filter_array_gpu!(ax_code_new, ax_code, idx, 0)

    @cusync begin
        ϕc = ϕc_new
        θc = θc_new
        μs = μs_new
        xx = xx_new
        yy = yy_new
        zz = zz_new
        wts = wts_new
        z_rot = z_rot_new
        ax_code = ax_code_new
    end

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
        allwavs = CUDA.zeros(precision, num_nonzero, 200)
        allints = CUDA.zeros(precision, num_nonzero, 200)
    end

    return GPUAllocs(λs_gpu, prof_gpu, ϕc, θc, μs, xx, yy, zz,
                     wts, z_rot, z_cbs, ax_code, dat_idx,
                     tloop_gpu, tloop_init, allwavs, allints)
end
