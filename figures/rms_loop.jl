# parallelized for loop
@everywhere function rms_loop(spec::SpecParams, disk::DiskParams, Nloop::T; top::Float64=NaN) where T<:Integer
    rvs = zeros(disk.Nt)
    rms0 = zeros(Nloop)
    for j in 1:Nloop
        rvs .= 0.0
        lambdas, outspec = synthesize_spectra(spec, disk, seed_rng=false, top=top)

        # extract velocities
        v_grid, ccf = calc_ccf(lambdas, outspec, spec, normalize=true)
        rvs, sigs = calc_rvs_from_ccf(v_grid, ccf)
        rms0[j] = calc_rms(rvs)
    end
    return mean(rms0), std(rms0)
end
