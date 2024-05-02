struct GPUAllocsEclipse{T1<:AF}
    λs::CuArray{T1,1}
    prof::CuArray{T1,1}
    flux::CuArray{T1,2}
end

function GPUAllocsEclipse(spec::SpecParams, disk::DiskParams; precision::DataType=Float64, verbose::Bool=true)
    # get dimensions for memory alloc
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    @cusync begin
        λs_gpu = CuArray{precision}(spec.lambdas)
        prof_gpu = CUDA.zeros(precision, Nλ)
        flux_gpu = CUDA.ones(precision, Nλ, Nt)
    end

    return GPUAllocs(λs_gpu, prof_gpu, flux_gpu)
end
