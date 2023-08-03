function generate_tloop!(tloop::AA{Int,2}, disk::DiskParams, soldata::SolarData{T}) where T<:AF
    # make sure dimensions are correct
    @assert size(tloop) == (length(disk.ϕc), length(disk.θc))

    # get the value of mu and ax codes
    disc_ax = parse_ax_string.(getindex.(keys(soldata.len),1))
    disc_mu = parse_mu_string.(getindex.(keys(soldata.len),2))

    # loop over grid
    xyz = zeros(3)
    for i in eachindex(disk.ϕc)
        for j in eachindex(disk.θc)
            # get cartesian coord
            xyz .= sphere_to_cart(disk.ρs, disk.ϕc[i], disk.θc[j])

            # calculate mu
            μc = calc_mu(xyz, disk.R_θ, disk.O⃗)

            # move to next iteration if patch element is not visible
            μc <= zero(T) && continue

            # get input data for place on disk
            key = get_key_for_pos(μc, xyz[1], xyz[3], disc_mu, disc_ax)
            len = soldata.len[key]

            # generate random index
            tloop[i,j] = floor(Int, rand() * len) + 1
        end
    end
    return nothing
end

"""
    synthesize_spectra(spec, disk; seed_rng=false, verbose=true, top=NaN)

Synthesize spectra given parameters in `spec` and `disk` instances.

# Arguments
- `spec::SpecParams`: SpecParams instance
- `disk::DiskParams`: DiskParams instance
"""
function synthesize_spectra(spec::SpecParams{T}, disk::DiskParams{T};
                            seed_rng::Bool=false, verbose::Bool=true,
                            use_gpu::Bool=false, precision::DataType=Float64,
                            skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # call appropriate simulation function on cpu or gpu
    if use_gpu
        return synth_gpu(spec, disk, seed_rng, verbose, precision, skip_times)
    else
        return synth_cpu(spec, disk, seed_rng, verbose, skip_times)
    end
end

function synth_cpu(spec::SpecParams{T}, disk::DiskParams{T}, seed_rng::Bool,
                   verbose::Bool, skip_times::BitVector) where T<:AF
    # parse out dimensions for memory allocation
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # allocate memory for time indices
    tloop = zeros(Int, length(disk.ϕc), length(disk.θc))
    tloop_init = similar(tloop)

    # allocate memory for spectra synthesis
    wsp = SynthWorkspace(ngrid=length(disk.ϕc))
    prof = ones(Nλ)
    outspec = ones(Nλ, Nt)
    outspec_temp = zeros(Nλ, Nt)


    # get number of calls to disk_sim needed
    templates = unique(spec.templates)

    # run the simulation (outspec modified in place)
    for (idx, file) in enumerate(templates)
        # re-seed the rng
        if seed_rng
            Random.seed!(42)
        end

        # get temporary specparams with lines for this run
        spec_temp = SpecParams(spec, file)

        # load in the appropriate input data
        if verbose
            println("\t>>> " * splitdir(file)[end])
        end
        soldata = SolarData(fname=file)

        # generate or copy tloop
        if (idx > 1) && in_same_group(templates[idx - 1], templates[idx])
            tloop .= tloop_init
        else
            generate_tloop!(tloop_init, disk, soldata)
            tloop .= tloop_init
        end

        # re-set array to 0s
        outspec_temp .= 0.0

        # run the simulation and multiply outspec by this spectrum
        disk_sim_3d(spec_temp, disk, soldata, wsp, prof, outspec_temp, tloop,
                    skip_times=skip_times, verbose=verbose)
        outspec .*= outspec_temp
    end
    return spec.lambdas, outspec
end

function synth_gpu(spec::SpecParams{T}, disk::DiskParams{T}, seed_rng::Bool,
                   verbose::Bool, precision::DataType, skip_times::BitVector) where T<:AF
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
    tloop_init = zeros(Int, length(disk.ϕc), length(disk.θc))
    outspec = ones(Nλ, Nt)

    # get number of calls to disk_sim needed
    templates = unique(spec.templates)

    # pre-allocate memory for gpu
    gpu_allocs = GPUAllocs(spec, disk, precision=precision)

    # run the simulation and return
    for (idx, file) in enumerate(templates)
        # re-seed the rng
        if seed_rng
            Random.seed!(42)
        end

        # get temporary specparams with lines for this run
        spec_temp = SpecParams(spec, file)

        # load in the appropriate input data
        if verbose
            println("\t>>> " * splitdir(file)[end])
        end
        soldata = SolarData(fname=file)

        # generate or copy tloop
        if (idx > 1) && in_same_group(templates[idx - 1], templates[idx])
            @cusync CUDA.copyto!(gpu_allocs.tloop, tloop_init)
        else
            generate_tloop!(tloop_init, disk, soldata)
            @cusync CUDA.copyto!(gpu_allocs.tloop, tloop_init)
        end

        # run the simulation and multiply outspec by this spectrum
        disk_sim_gpu(spec_temp, disk, soldata, gpu_allocs, outspec, verbose=verbose,
                     skip_times=skip_times, precision=precision)
    end
    return spec.lambdas, outspec
end
