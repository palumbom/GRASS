using Pkg; Pkg.activate(".")
using JLD2
using CUDA
using GRASS
using Printf
using FileIO
using Statistics
using BenchmarkTools

# plotting imports
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")

# parse args + get directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

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
    max_cpu = minimum([12, length(lines)])
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
    b_gpu_avg = dropdims(mean(b_gpu, dims=2),dims=2)
    b_gpu_std = dropdims(std(b_gpu, dims=2),dims=2)

    # compute speedup
    speedup = b_cpu[1:max_cpu]./b_gpu_avg[1:max_cpu]

    println(">>> Max GPU benchmark = " * string(maximum(b_gpu_avg)))

    # plotting function (use globals who cares i don't)
    function plot_scaling(;logscale=true)
        # create plotting objects
        fig, ax1 = plt.subplots()
        ax2 = ax1.twiny()

        # log scale it
        if logscale
            ax1.set_yscale("symlog")
            ax2.set_yscale("symlog")
            scale = "logscale"
        else
            scale = "linscale"
        end

        # plot on ax1
        ms = 7.5
        ax1.plot(n_res[1:max_cpu], b_cpu[1:max_cpu], marker="o", ms=ms, c="k", label=L"{\rm CPU}")
        ax1.plot(n_res, b_gpu_avg, marker="s", ms=ms, c="tab:blue", label=L"{\rm GPU}")

        # plot on twin axis
        ax2.plot(n_lam[1:max_cpu], b_cpu[1:max_cpu], marker="o", ms=ms, c="k")
        ax2.plot(n_lam, b_gpu_avg, marker="s", ms=ms, c="tab:blue")
        ax2.grid(false)

        # axis label stuff
        ax1.set_xlabel(L"{\rm \#\ of\ res.\ elements}")
        ax1.set_ylabel(L"{\rm Synthesis\ time\ (s)}")
        ax2.set_xlabel(L"{\rm Width\ of\ spectrum\ (\AA)}")
        ax1.legend()
        plt.show()
        fig.savefig(plotdir * "scaling_bench_" * scale * ".pdf")
        plt.clf(); plt.close()
        return nothing
    end

    # plot it
    plot_scaling(logscale=true)
    plot_scaling(logscale=false)
end
