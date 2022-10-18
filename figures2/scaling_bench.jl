using Pkg; Pkg.activate(".")
using Printf
using Profile
using BenchmarkTools
using CUDA
using GRASS

using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
# mpl.style.use("my.mplstyle")

# define some functions
grassdir, plotdir, datadir = check_plot_dirs()

@assert CUDA.functional()

# benchmark on cpu
function benchmark_cpu(lines::AbstractArray{Float64,1}, depths::AbstractArray{Float64,1})
    templates = repeat(["FeI_5434"], length(lines))
    resolution = 700000.0
    spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution)
    disk = DiskParams(N=132, Nt=50)

    # synthesize it
    lambdas1, outspec1 = synthesize_spectra(spec, disk, verbose=false)
    return nothing
end

# benchmark on gpu
function benchmark_gpu(lines::AbstractArray{Float64,1}, depths::AbstractArray{Float64,1})
    templates = repeat(["FeI_5434"], length(lines))
    resolution = 700000.0
    spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution)
    disk = DiskParams(N=132, Nt=50)

    # synthesize it
    lambdas1, outspec1 = synthesize_spectra(spec, disk, verbose=false, use_gpu=true)
    return nothing
end

# line parameters
nlines = 24
lines = range(5434.5, step=5.0, length=nlines)
depths = repeat([0.75], length(lines))

# initial for loop to calculate number of reoslution elements per spectrum
n_res = similar(lines)
n_lam = similar(lines)
for i in eachindex(lines)
    # calculate coverage
    lines1 = lines[1:i]
    coverage = (minimum(lines1) - 0.75, maximum(lines1) + 0.75)

    # generate Delta ln lambda
    Δlnλ = (1.0 / 700000.0)
    lnλs = range(log(coverage[1]), log(coverage[2]), step=Δlnλ)
    lambdas = exp.(lnλs)
    n_res[i] = length(lambdas)
    n_lam[i] = lambdas[end] - lambdas[1]
end

# benchmarking wrapper function
function bmark_everything(b_cpu, b_gpu, lines, depths; max_cpu=8)
    # initial run to precompile things
    spec1 = SpecParams(lines=[5434.5], depths=[0.75], templates=["FeI_5434"], resolution=7e5)
    disk1 = DiskParams(N=132, Nt=10)
    lambdas1, outspec1 = synthesize_spectra(spec1, disk1, seed_rng=false, use_gpu=false)
    lambdas1, outspec1 = synthesize_spectra(spec1, disk1, seed_rng=false, use_gpu=true)

    for i in eachindex(lines)
        @printf("Benchmarking %s of %s\n", i, length(lines))

        lines_i = lines[1:i]
        depths_i = depths[1:i]

        if i <= max_cpu
            @printf("Performing CPU bench (N=132, Nt=50): ")
            Profile.clear_malloc_data()
            b_cpu[i] = @belapsed benchmark_cpu($lines_i, $depths_i)
            println()
        end

        @printf("Performing GPU bench (N=132, Nt=50): ")
        Profile.clear_malloc_data()
        b_gpu[i] = @belapsed benchmark_gpu($lines_i, $depths_i)
        println()
    end
    return nothing
end

# allocate memory for benchmark results and run it
num = 8
max_cpu = minimum([num, length(lines)])
b_cpu = similar(lines)
b_gpu = similar(lines)
bmark_everything(b_cpu, b_gpu, lines, depths, max_cpu=max_cpu)

# function for plotting since it dies in global??
function plot_it()
    fig, ax1 = plt.subplots()

    # plot on ax1
    ax1.plot(n_res[1:max_cpu], b_cpu[1:max_cpu], c="k", label="CPU")
    ax1.scatter(n_res[1:max_cpu], b_cpu[1:max_cpu], c="k")
    ax1.plot(n_res, b_gpu, c="tab:blue", label="GPU")
    ax1.scatter(n_res, b_gpu, c="tab:blue")

    # plot on twin axis
    ax2 = ax1.twiny()
    ax2.plot(n_lam[1:max_cpu], b_cpu[1:max_cpu], c="k")
    ax2.scatter(n_lam[1:max_cpu], b_cpu[1:max_cpu], c="k")
    ax2.plot(n_lam, b_gpu, c="tab:blue")
    ax2.scatter(n_lam, b_gpu, c="tab:blue")

    ax1.set_xlabel("Number of spectral resolution elements")
    ax1.set_ylabel("Simulation time (s)")
    ax2.set_xlabel("Width of spectrum (A)")
    ax1.legend()
    fig.savefig(plotdir * "scaling_bench.pdf")
    return nothing
end

plot_it()

