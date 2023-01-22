using Pkg; Pkg.activate(".")
using Printf
using Profile
using BenchmarkTools
using CUDA
using GRASS
using JLD2
using FileIO
using Statistics

# parse args + get directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

function precompile_with_spectra_test()
    N = 132
    Nt = 50
    lines = [5434.5]
    depths = [0.5]
    geffs = [0.0]
    templates = ["FeI_5434"]
    variability = repeat([true], length(lines))
    resolution = 7e5

    disk = DiskParams(N=N, Nt=Nt)
    spec = SpecParams(lines=lines, depths=depths, variability=variability,
                       geffs=geffs, templates=templates, resolution=resolution)

    lambdas1, outspec1 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=false)
    lambdas2, outspec2 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=true)
    return nothing
end

function benchmark_cpu(spec::SpecParams, disk::DiskParams)
    lambdas1, outspec1 = synthesize_spectra(spec, disk, verbose=false, seed_rng=true)
    return nothing
end

function benchmark_gpu(spec::SpecParams, disk::DiskParams)
    lambdas1, outspec1 = synthesize_spectra(spec, disk, verbose=false, seed_rng=true, use_gpu=true)
    return nothing
end

# benchmarking wrapper function
function bmark_everything(b_cpu, b_gpu, lines, depths; max_cpu=8)
    for i in eachindex(lines)
        @printf(">>> Benchmarking %s of %s\n", i, length(lines))

        # get number of gpu samplings to do
        n_gpu_loops = size(b_gpu, 2)

        # get lines and depths
        lines_i = lines[1:i]
        depths_i = depths[1:i]
        templates_i = repeat(["FeI_5434"], length(lines_i))
        resolution_i = 7e5

        # create spec, disk instances
        spec = SpecParams(lines=lines_i, depths=depths_i, templates=templates_i, resolution=resolution_i)
        disk = DiskParams(N=132, Nt=50)

        # CPU bench loop
        if i <= max_cpu
            @printf("\t>>> Performing CPU bench (N=132, Nt=50): ")
            Profile.clear_malloc_data()
            b_cpu[i] = @belapsed benchmark_cpu($spec, $disk)
            println()
        end

        # GPU bench loop
        @printf("\t>>> Performing GPU bench (N=132, Nt=50): ")
        for j in 1:n_gpu_loops
            Profile.clear_malloc_data()
            CUDA.synchronize()
            b_gpu[i,j] = @belapsed benchmark_gpu($spec, $disk)
            CUDA.synchronize()
        end
        println()
    end
    return nothing
end

function main()
    # make sure its precompiled
    println(">>> Doing precompilation synthesis...")
    precompile_with_spectra_test()

    # line parameters
    nlines = 24
    lines = range(5434.5, step=5.0, length=nlines)
    depths = repeat([0.75], length(lines))

    # calculate number of resolution elements per spectrum
    n_res = similar(lines)
    n_lam = similar(lines)
    for i in eachindex(lines)
        # calculate wavelength coverage
        lines1 = lines[1:i]
        coverage = (minimum(lines1) - 0.75, maximum(lines1) + 0.75)

        # generate Delta ln lambda
        Δlnλ = (1.0 / 7e5)
        lnλs = range(log(coverage[1]), log(coverage[2]), step=Δlnλ)
        lambdas = exp.(lnλs)
        n_res[i] = length(lambdas)
        n_lam[i] = last(lambdas) - first(lambdas)
    end

    # allocate memory for benchmark results and run it
    n_gpu_loops = 16
    max_cpu = minimum([8, length(lines)])
    b_cpu = similar(lines)
    b_gpu = zeros(length(lines), n_gpu_loops)
    bmark_everything(b_cpu, b_gpu, lines, depths, max_cpu=max_cpu)

    # write info to disk
    outfile = datadir * "scaling_benchmark.jld2"
    save(outfile,
         "max_cpu", max_cpu,
         "nlines", nlines,
         "n_res", n_res,
         "n_lam", n_lam,
         "b_cpu", b_cpu,
         "b_gpu", b_gpu)
    return nothing
end

if run
    main()
end

if plot
    # plotting imports
    using LaTeXStrings
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

    # read in the data
    file = datadir * "scaling_benchmark.jld2"
    d = load(file)
    max_cpu = d["max_cpu"]
    nlines = d["nlines"]
    n_res = d["n_res"]
    n_lam = d["n_lam"]
    b_cpu = d["b_cpu"]
    b_gpu = d["b_gpu"]

    # get mean gpu benchmark
    b_gpu_avg = mean(b_gpu, dims=2)
    b_gpu_std = std(b_gpu, dims=2)

    # create plotting objects
    fig, ax1 = plt.subplots()
    ax2 = ax1.twiny()
    # ax3 = ax1.twinx()

    # log scale it
    # ax1.set_yscale("symlog")
    # ax2.set_yscale("symlog")
    # ax2.set_yscale("symlog")


    # plot on ax1
    ms = 5.0
    ax1.plot(n_res[1:max_cpu], b_cpu[1:max_cpu], marker="o", ms=ms, c="k", label=L"{\rm CPU}")
    ax1.plot(n_res, b_gpu_avg, marker="s", ms=ms, c="tab:blue", label=L"{\rm GPU}")

    # plot on twin axis
    ax2.plot(n_lam[1:max_cpu], b_cpu[1:max_cpu], marker="o", ms=ms, c="k")
    ax2.plot(n_lam, b_gpu_avg, marker="s", ms=ms, c="tab:blue")
    ax2.grid(false)

    # plot on twin axis
    # ax3.plot(n_lam[1:max_cpu], b_cpu[1:max_cpu], marker="o", ms=ms, c="k")
    # ax3.plot(n_lam, b_gpu_avg, marker="s", ms=ms, c="tab:blue")
    # ax3.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax3.get_yaxis().set_minor_formatter(mpl.ticker.NullFormatter())
    # ax3.grid(false)

    ax1.set_xlabel(L"{\rm \#\ of\ res.\ elements}")
    ax1.set_ylabel(L"{\rm Synthesis\ time\ (s)}")
    ax2.set_xlabel(L"{\rm Width\ of\ spectrum\ (\AA)}")
    ax1.legend()
    fig.savefig(plotdir * "scaling_bench.pdf")
    plt.clf(); plt.close()
end
