"""
    synthesize_spectra_eclipse(
        spec,
        disk,
        wavelength,
        LD_type,
        obs_long,
        obs_lat,
        alt,
        time_stamps,
        ext_coeff;
        ext_toggle=false,
        seed_rng=false,
        verbose=true,
        use_gpu=false,
        precision=Float64,
        skip_times=falses(disk.Nt),
    )

Synthesize a time series of solar spectra as seen from a ground-based observer
during a solar eclipse. For each epoch in `time_stamps`, the Sun–Moon–observer
geometry is computed from SPICE ephemerides, the lunar disk is projected onto the
solar grid to mask occulted cells, and the remaining surface is integrated to
produce a spectrum. Returns a tuple `(λs, flux)` where `λs::Vector` is the
wavelength grid (`spec.lambdas`) and `flux::Matrix` has size `(length(λs), disk.Nt)`,
one column per epoch.

Required SPICE kernels are fetched and furnished on the first call.

# Arguments
- `spec::SpecParams`: spectral synthesis parameters (line list, templates, wavelength grid).
- `disk::DiskParamsEclipse`: eclipse disk parameters; see [`DiskParamsEclipse`](@ref).
- `wavelength::Vector{Float64}`: line-center wavelength(s) used to look up limb-darkening coefficients.
- `LD_type::String`: limb-darkening law. One of `"SSD_quadratic"`, `"SSD_4parameter"`, `"300_quadratic"`, `"300_4parameter"`, `"HD_quadratic"`, `"HD_4parameter"` (dataset `SSD`/`300`/`HD` × quadratic or 4-parameter law).
- `obs_long::T`, `obs_lat::T`: observer geodetic longitude and latitude (degrees).
- `alt::T`: observer altitude (km).
- `time_stamps::Vector{String}`: UTC timestamp strings (SPICE-parseable), one per epoch. Length must equal `disk.Nt`.
- `ext_coeff`: atmospheric extinction coefficient (used only when `ext_toggle=true`).

# Keyword Arguments
- `ext_toggle::Bool=false`: include airmass-dependent atmospheric extinction.
- `seed_rng::Bool=false`: re-seed the RNG with a fixed seed per template for reproducible granulation.
- `verbose::Bool=true`: print progress messages for template loading and simulation.
- `use_gpu::Bool=false`: run the GPU implementation (requires a functional CUDA device).
- `precision::DataType=Float64`: GPU precision (`Float32` or `Float64`); see [Caveats](@ref "Caveats").
- `skip_times::BitVector=falses(disk.Nt)`: epochs to skip in the simulation loop.
"""
function synthesize_spectra_eclipse(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, wavelength::Vector{Float64}, LD_type::String, 
                                    obs_long::T, obs_lat::T, alt::T, time_stamps::Vector{String},                              
                                    ext_coeff; ext_toggle::Bool=false, seed_rng::Bool=false, verbose::Bool=true,
                                    use_gpu::Bool=false, precision::DataType=Float64,
                                    skip_times::BitVector=falses(disk.Nt)) where T<:AF
    GRASS.Eclipse.get_kernels()
    
    # call appropriate simulation function on cpu or gpu
    if use_gpu
        return synth_Eclipse_gpu(spec, disk, verbose, precision, skip_times, LD_type,
                                    obs_long, obs_lat, alt, time_stamps, wavelength, ext_coeff, ext_toggle, false)
    else
        return synth_Eclipse_cpu(spec, disk, seed_rng, verbose, skip_times, LD_type, wavelength, 
                                    time_stamps, obs_long, obs_lat, alt, ext_coeff, ext_toggle)
    end
end

function synth_Eclipse_cpu(spec::SpecParams{T}, disk::DiskParamsEclipse{T}, seed_rng::Bool,
                            verbose::Bool, skip_times::BitVector, LD_type::String, wavelength::Vector{Float64}, 
                            time_stamps::Vector{String}, obs_long::T, obs_lat::T, alt::T, ext_coeff, ext_toggle::Bool) where T<:AF

    # parse out dimensions for memory allocation
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # allocate memory for synthsis
    prof = ones(Nλ)
    flux = ones(Nλ, Nt)

    # pre-allocate memory and pre-compute geometric quantities
    wsp = SynthWorkspaceEclipse(disk, Int(length(wavelength)), Nt, verbose=verbose)

    # allocate memory for time indices
    tloop = zeros(Int, size(disk.θc))
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
        soldata = SolarData(fname=file)

        # re-seed the rng
        if seed_rng
            Random.seed!(42)
        end

        # run the simulation and multiply flux by this spectrum
        GRASS.Eclipse.disk_sim_eclipse(spec_temp, disk, soldata, wsp, prof, flux, tloop, tloop_init, templates, idx, LD_type, wavelength, 
                        time_stamps, obs_long, obs_lat, alt, ext_coeff, ext_toggle, skip_times=skip_times)
    end
    return spec.lambdas, flux
end

function synth_Eclipse_gpu(spec::SpecParams{T}, disk::DiskParamsEclipse{T},
                            verbose::Bool, precision::DataType, skip_times::BitVector, LD_type::String, 
                            obs_long::T, obs_lat::T, alt::T, time_stamps::Vector{String}, 
                            wavelength, ext_coeff, ext_toggle::Bool, spot_toggle::Bool) where T<:AF
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
        soldata_cpu = SolarData(fname=file)
        soldata = GPUSolarData(soldata_cpu, precision=precision)

        # run the simulation and multiply flux by this spectrum
        GRASS.Eclipse.disk_sim_eclipse_gpu(spec_temp, disk, soldata, gpu_allocs, flux, 
                              obs_long, obs_lat, alt, time_stamps, wavelength, 
                              ext_coeff, ext_toggle, spot_toggle, LD_type, skip_times=skip_times)
    end
    return spec.lambdas, flux
end

function synth_Eclipse_gpu(spec::SpecParams{T}, disk::DiskParamsEclipse{T},
                            verbose::Bool, precision::DataType, skip_times::BitVector, 
                            obs_long::T, obs_lat::T, alt::T, time_stamps::Vector{String}, 
                            wavelength, ext_coeff, CB1, CB2, CB3) where T<:AF
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
        soldata_cpu = SolarData(fname=file)
        soldata = GPUSolarData(soldata_cpu, precision=precision)

        # run the simulation and multiply flux by this spectrum
        GRASS.Eclipse.disk_sim_eclipse_gpu(spec_temp, disk, soldata, gpu_allocs, flux, 
                              obs_long, obs_lat, alt, time_stamps, wavelength, 
                              ext_coeff, CB1, CB2, CB3, skip_times=skip_times)
    end
    return spec.lambdas, flux
end