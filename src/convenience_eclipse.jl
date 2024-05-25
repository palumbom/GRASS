"""
    synthesize_spectra(spec, disk; seed_rng=false, verbose=true, top=NaN)

Synthesize spectra given parameters in `spec` and `disk` instances.

# Arguments
- `spec::SpecParams`: SpecParams instance
- `disk::DiskParams`: DiskParams instance
"""
function synthesize_spectra_eclipse(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, obs_long, obs_lat, alt, wavelength, time_stamps;
                            seed_rng::Bool=false, verbose::Bool=true,
                            use_gpu::Bool=false, precision::DataType=Float64, fixed_bisector::Bool=false,
                            skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # call appropriate simulation function on cpu or gpu
    if use_gpu
        return synth_Eclipse_gpu(spec, disk, seed_rng, verbose, precision, skip_times, obs_long, obs_lat, alt, time_stamps, wavelength, fixed_bisector)
    else
        return synth_Eclipse_cpu(spec, disk, seed_rng, verbose, skip_times, obs_long, obs_lat, alt, time_stamps, wavelength, fixed_bisector)
    end
end

function synth_Eclipse_cpu(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, seed_rng::Bool,
                   verbose::Bool, skip_times::BitVector, obs_long, obs_lat, alt, time_stamps, wavelength, fixed_bisector::Bool) where T<:AF

    # parse out dimensions for memory allocation
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # allocate memory for synthsis
    prof = ones(Nλ)
    flux = ones(Nλ, Nt)

    # pre-allocate memory and pre-compute geometric quantities
    wsp = SynthWorkspaceEclipse(disk, Int(length(wavelength)), verbose=verbose)
    mem = GeoWorkspaceEclipse(disk, Int(length(wavelength)))

    # allocate memory for time indices
    tloop = zeros(Int, size(wsp.μs))
    tloop_init = zeros(Int, size(tloop))

    # get number of calls to disk_sim needed
    templates = unique(spec.templates)

    # run the simulation (flux modified in place)
    for (idx, file) in enumerate(templates)
        # get temporary specparams with lines for this run
        spec_temp = SpecParams(spec, file)

        # load in the appropriate input data
        if verbose
            println("\t>>> Template: " * splitdir(file)[end])
        end
        soldata = SolarData(fname=file, fixed_bisector=fixed_bisector)

        # re-seed the rng
        if seed_rng
            Random.seed!(42)
        end

        # run the simulation and multiply flux by this spectrum
        disk_sim_eclipse(spec_temp, disk, soldata, wsp, mem, prof, flux, tloop, tloop_init, templates, idx, obs_long, obs_lat, alt, time_stamps, wavelength,
                 fixed_bisector, skip_times=skip_times, verbose=verbose)
    end
    return spec.lambdas, flux
end

function synth_Eclipse_gpu(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, seed_rng::Bool,
                   verbose::Bool, precision::DataType, skip_times::BitVector, obs_long, obs_lat, alt, time_stamps, wavelength, fixed_bisector::Bool) where T<:AF
    # make sure there is actually a GPU to use
    @assert CUDA.functional()

    # warn user if precision is single
    if precision <: Float32
       @warn "Single-precision GPU implementation produces large flux and velocity errors!"
    end

    # parse out dimensions for memory allocation
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # allocate memory
    flux = ones(Nλ, Nt)

    # get number of calls to disk_sim needed
    templates = unique(spec.templates)

    # allocate memory needed for rossiter computations
    gpu_allocs = GPUAllocsEclipse(spec, disk, Int(length(wavelength)), precision=precision, verbose=verbose)

    # run the simulation and return
    for (idx, file) in enumerate(templates)
        # get temporary specparams with lines for this run
        spec_temp = SpecParams(spec, file)

        # load in the appropriate input data
        if verbose
            println("\t>>> Template: " * splitdir(file)[end])
        end
        soldata_cpu = SolarData(fname=file, fixed_bisector=fixed_bisector)
        soldata = GPUSolarData(soldata_cpu, precision=precision)

        # run the simulation and multiply flux by this spectrum
        disk_sim_eclipse_gpu(spec_temp, disk, soldata, gpu_allocs,
                              flux, templates, idx, 
                              obs_long, obs_lat, alt, time_stamps, wavelength, verbose=verbose,
                              skip_times=skip_times)
    end
    return spec.lambdas, flux
end