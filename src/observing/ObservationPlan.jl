const def_exp_tim = 15

struct ObservationPlan{T1<:Integer, T2<:AbstractFloat}
    N_obs::T1
    N_nights::T1
    obs_per_night::T1
    time_per_obs::T2
    dead_time::T2
end

function ObservationPlan(;N_obs::T1=0, obs_per_night::T1=0, time_per_obs::T2=0.0,
                         dead_time::T2=0.0) where {T1<:Integer, T2<:AF}
    # assertions
    @assert N_obs > 0
    @assert obs_per_night > 0
    @assert time_per_obs >= def_exp_tim
    @assert dead_time >= 0.0

    # TODO only works if exposure times are 0 or factor of 15
    @assert iszero(time_per_obs % 15)
    @assert iszero(dead_time % 15)

    # get number of nights
    N_nights = Int(N_obs ÷ obs_per_night) + one(T1) * !iszero(N_obs % obs_per_night)
    return ObservationPlan(N_obs, N_nights, obs_per_night, time_per_obs, dead_time)
end

function get_nt(obs::ObservationPlan)
    nt_per_exposure = Int(obs.time_per_obs ÷ def_exp_tim) + 1 * !iszero(obs.time_per_obs % def_exp_tim)
    nt_per_deadtime = Int(obs.dead_time ÷ def_exp_tim) + 1 * !iszero(obs.dead_time % def_exp_tim)
    return nt_per_exposure, nt_per_deadtime
end

function get_skip_times(nt_per_night, nt_per_exposure, nt_per_deadtime)
    # get relevant numbers for iteration
    nt_per_iter = nt_per_exposure + nt_per_deadtime

    # idk fuck around
    iter_skips = BitArray{1}(undef, nt_per_iter)
    iter_skips[1:nt_per_exposure] .= false
    iter_skips[nt_per_exposure+1:end] .= true

    # allocate memory
    return repeat(iter_skips, nt_per_night ÷ nt_per_iter)
end

function simulate_observations(obs::ObservationPlan, spec::SpecParams;
                               N::Int=128, snr::T=Inf, new_res::T=NaN,
                               use_gpu::Bool=false) where T<:AF
    # assertions
    @assert (isnan(new_res) | (new_res < 7e5))

    # get number of time steps needed per exposure
    nt_per_exposure, nt_per_deadtime = get_nt(obs)

    # get number of time steps needed per night
    nt_per_iter = nt_per_exposure + nt_per_deadtime
    nt_per_night = obs.obs_per_night * nt_per_iter

    # get BitArray of epochs to skip
    skip_times = get_skip_times(nt_per_night, nt_per_exposure, nt_per_deadtime)

    # set up disk params instance
    disk = DiskParams(N=N, Nt=nt_per_night)

    # synthesize the spectra at default resolution
    Nλ = length(spec.lambdas)
    # outspec = zeros(Nλ, disk.Nt)

    wavs, flux = synthesize_spectra(spec, disk, verbose=false, use_gpu=use_gpu,
                                    skip_times=skip_times)

    # if use_gpu
    #     disk_sim_gpu(spec, disk, outspec, skip_times=skip_times)
    # else
    #     disk_sim(spec, disk, prof, outspec, skip_times=skip_times)
    #     prof = ones(Nλ)
    # end


    # allocate memory for binned spectra
    flux_binned = zeros(Nλ, obs.N_obs)
    for i in 1:obs.N_obs
        inds = (i-1) * nt_per_iter + 1: i * nt_per_iter
        flux_binned[:, i] = sum(flux[:, inds], dims=2) ./ sum(.!skip_times[inds])
    end

    # degrade the resolution
    if !isnan(new_res)
        # set two pixels per res element
        oversampling = 2.0

        # get size of output to convolve from initial convolution
        wavs_to_deg = view(wavs, :, 1)
        flux_to_deg = view(flux_binned, :, 1)
        wavs_degd, flux_degd = GRASS.convolve_gauss(wavs_to_deg,
                                                    flux_to_deg,
                                                    new_res=new_res,
                                                    oversampling=oversampling)

        # allocate memory for degraded spectra
        # wavs_degd = zeros(size(wavs_degd, 1), size(flux_binned, 2))
        flux_degd = zeros(size(wavs_degd, 1), size(flux_binned, 2))

        # loop over epochs and do convolution
        for j in 1:size(flux_binned, 2)
            flux_to_deg = view(flux_binned, :, j)
            wavs_temp, flux_temp = GRASS.convolve_gauss(wavs_to_deg,
                                                        flux_to_deg,
                                                        new_res=new_res,
                                                        oversampling=oversampling)

            # copy to array
            # wavs_degd[:, j] .= wavs_temp
            flux_degd[:, j] .= flux_temp
        end

        wavs = wavs_degd
        flux_binned = flux_degd
    end

    # add noise to specified snr per res element
    flux_binned = add_noise(flux_binned, snr)

    # calculate ccf and velocities
    v_grid, ccf1 = calc_ccf(wavs, flux_binned, spec, normalize=true)
    rvs, sigs = calc_rvs_from_ccf(v_grid, ccf1)
    return rvs, sigs
end
