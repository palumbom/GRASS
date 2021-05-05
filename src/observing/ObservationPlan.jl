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
    @assert N_obs > 0
    @assert obs_per_night > 0
    @assert time_per_obs >= def_exp_tim
    @assert dead_time >= 0.0

    # TODO this only works if dead time is 0 or factor of 15
    # @assert iszero(mod(dead_time, def_exp_tim))

    # warn if exposure time is not factor of 15
    if (!iszero(time_per_obs % def_exp_tim) & !iszero((time_per_obs + dead_time) % def_exp_tim))
        println(">>> Warning: Implementation for exposure times that aren't a factor of 15 is funky! \n")
    end

    # get number of nights
    N_nights = Int(N_obs ÷ obs_per_night) + one(T1) * !iszero(N_obs % obs_per_night)
    return ObservationPlan(N_obs, N_nights, obs_per_night, time_per_obs, dead_time)
end

function get_nt_from_time(obs::ObservationPlan, t::Symbol)
    @assert t in fieldnames(ObservationPlan)
    return Int(getfield(obs, t) ÷ def_exp_tim) + 1 * !iszero(getfield(obs, t) % def_exp_tim)
end

function get_nt_from_time(t::T) where T<:Real
    return Int(t ÷ def_exp_tim) + 1 * !iszero(t % def_exp_tim)
end

function simulate_observations(obs::ObservationPlan, spec::SpecParams;
                               N::Int=128, snr::T=Inf, new_res::T=NaN) where T<:AF
    @assert (isnan(new_res) | (new_res < 7e5))

    # set boolean for long dead time
    if obs.dead_time <= (60.0 * 20.0) # less than 20 minutes
        long_dead = false
    else
        long_dead = true
    end

    # get number of time steps needed per exposure
    nt_per_exposure = get_nt_from_time(obs, :time_per_obs)
    nt_per_deadtime = get_nt_from_time(obs, :dead_time)

    # get number of time steps needed per night
    nt_per_iter = nt_per_exposure + nt_per_deadtime
    nt_per_night = obs.obs_per_night * (nt_per_exposure + nt_per_deadtime)

    # get array of weights
    weights = zeros(nt_per_night)
    nt = 0
    dead = false
    for i in eachindex(weights)
        # initial increment of nt
        nt += 1

        # set weight for exposure
        if ((nt <= nt_per_exposure) & !dead)

            # TODO this line is messed up; figure it out
            num = (obs.time_per_obs - def_exp_tim * (nt - 1))
            if num > def_exp_tim
                num = def_exp_tim
            end
            weights[i] = num / def_exp_tim
            if nt == nt_per_exposure
                nt = 0
                if !iszero(nt_per_deadtime)
                    dead = true
                end
            end
            continue
        end

        # set weight to 0 for dead time
        # TODO: broken for fraction cases
        if (!iszero(nt_per_deadtime) & (nt <= nt_per_deadtime) & dead)
            weights[i] = NaN
            if nt == nt_per_deadtime
                nt = 0
                dead = false
            end
        end
    end

    # allocate memory for binning spectra and for velocities
    vels = zeros(obs.N_obs)
    flux_bin = zeros(length(spec.lambdas), obs.obs_per_night)

    # loop over the nights
    for i in 0:(obs.N_nights-1)
        # simulate the spectra if dead time is small
        if !long_dead
            # get disk params from num time and synthesize spectra
            disk = DiskParams(N=N, Nt=nt_per_night)
            wavs, flux = synthesize_spectra(spec, disk, top=NaN)
        end

        # loop over observations within night
        for j in 0:(obs.obs_per_night-1)
            # simulate the spectra if dead time is large
            if long_dead
                disk = DiskParams(N=N, Nt=nt_per_exposure)
                wavs, flux = synthesize_spectra(spec, disk, top=NaN)
                itr = 1:nt_per_exposure
            else
                itr = (j*nt_per_iter+1):(j*nt_per_iter + nt_per_iter)
            end

            # multiply by the weights
            for k in itr
                flux[:, k] .*= weights[k]
            end

            # do the binning on non-NaN values
            flux_to_bin = strip_nans_by_column(flux[:, itr])

            # TODO: nanmath package?
            flux_bin[:, j+1] = sum(flux_to_bin, dims=2)/sum(filter(!isnan, weights[itr]))
        end

        # add noise to specificed snr
        # TODO: noise per pixel vs resolution element?
        flux_bin = add_noise(flux_bin, snr)

        # degrade the resolution
        if !isnan(new_res)
            flux_bin = convolve_gauss(wavs, flux_bin, new_res=new_res)
        end

        # calculate velocities
        v_grid, ccf1 = calc_ccf(wavs, flux_bin, spec, normalize=true)
        rvs, sigs = calc_rvs_from_ccf(v_grid, ccf1)

        # copy the rvs to the output array
        itr = (i*obs.obs_per_night+1):(i*obs.obs_per_night + obs.obs_per_night)
        vels[itr] = rvs
    end
    return vels
end
