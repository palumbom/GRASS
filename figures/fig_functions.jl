using Distributed

# parallelized for loop
@everywhere function spec_loop(spec::SpecParams, disk::DiskParams, Nloop::T; top::Float64=NaN) where T<:Integer
    rms0 = zeros(Nloop)
    avg0 = zeros(Nloop)
    for j in 1:Nloop
        # synthesize spectra
        lambdas, outspec = synthesize_spectra(spec, disk, seed_rng=false, top=top)

        # extract velocities
        v_grid, ccf = calc_ccf(lambdas, outspec, spec, normalize=true)
        rvs, sigs = calc_rvs_from_ccf(v_grid, ccf)

        # compute rms and mean
        rms0[j] = calc_rms(rvs)
        avg0[j] = mean(rvs)
    end
    return mean(avg0), std(avg0), mean(rms0), std(rms0)
end

# small function to check and create file structure
function check_plot_dirs()
    # directories
    dirs = [homedir() * "/grass_output/",
            homedir() * "/grass_output/plots/",
            homedir() * "/grass_output/data/"]

    # create dirs if they dont exist, and return dir names
    for dir in dirs
        if !isdir(dir)
            mkdir(dir)
        end
    end
    return dirs
end

# small function to parse command line arguments
function parse_args(args)
    # chain of if/else to set booleans
    if isempty(args)
        run = true
        plot = true
    else
        if "run" in args
            run = true
        else
            run = false
        end

        if "plot" in args
            plot = true
        else
            plot = false
        end
    end

    # print statements for clarity
    if (run & plot)
        println(">>> Running and plotting " * string(PROGRAM_FILE) * "...")
    elseif run
        println(">>> Running " * string(PROGRAM_FILE) * "...")
    elseif plot
        println(">>> Plotting " * string(PROGRAM_FILE) * "...")
    else
        println(">>> Only defining functions in " * string(PROGRAM_FILE) * "...")
    end
    return run, plot
end
