# parallelized for loop
function spec_loop(spec::SpecParams, disk::DiskParams, Nloop::T; use_gpu::Bool=false) where T<:Integer
    rms0 = zeros(Nloop)
    avg0 = zeros(Nloop)
    for j in 1:Nloop
        # synthesize spectra
        lambdas, outspec = synthesize_spectra(spec, disk, seed_rng=false, use_gpu=use_gpu, verbose=false)
        CUDA.memory_status()
        CUDA.reclaim()
        println()

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
function check_plot_dirs(;topdir=nothing)
    # set the directory you want the outpout to be in
    if isnothing(topdir)
        topdir = homedir()
    end

    # directories
    dirs = [topdir * "/grass_output/",
            topdir * "/grass_output/plots/",
            topdir * "/grass_output/data/"]

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

