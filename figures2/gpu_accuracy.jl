using Pkg; Pkg.activate(".")
using JLD2
using CUDA
using GRASS
using Printf
using FileIO
using Revise
using Statistics
using EchelleCCFs
using BenchmarkTools

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

# parse args + get directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# function to create spectra
function generate_spectrum(disk::DiskParams, spec::SpecParams; use_gpu::Bool=false)
    use_gpu && @assert CUDA.functional()
    wavs, flux = synthesize_spectra(spec, disk, seed_rng=true, verbose=true, use_gpu=use_gpu)
    return wavs, flux
end

# set up paramaters for spectrum
lines = [5500.0, 5501.5]
depths = [0.8, 0.65]
templates = ["FeI_5250.6", "FeI_5434"]
resolution = 7e5
buffer = 1.5

# create composite types
disk = DiskParams(N=132, Nt=5)
spec = SpecParams(lines=lines, depths=depths, templates=templates,
                  resolution=resolution, buffer=buffer)

# get the spectra
wavs_cpu, flux_cpu = generate_spectrum(disk, spec, use_gpu=false)
wavs_gpu, flux_gpu = generate_spectrum(disk, spec, use_gpu=true)

# compute means
flux_cpu_mean = dropdims(mean(flux_cpu, dims=2), dims=2)
flux_gpu_mean = dropdims(mean(flux_gpu, dims=2), dims=2)

# get residuals
resids = flux_gpu_mean .- flux_cpu_mean

# set up plot
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8,6), sharex=true)
fig.subplots_adjust(hspace=-0.05)
# ax2.ticklabel_format(style="scientific", axis="y")

ax1.plot(wavs_cpu, flux_cpu_mean, ls="--", c=colors[1], label=L"{\rm CPU}")
ax1.plot(wavs_gpu, flux_gpu_mean, ls="-.", c=colors[2], label=L"{\rm GPU}")
ax2.scatter(wavs_cpu, resids, s=2, c="k")

ax1.set_ylabel(L"{\rm Normalized\ Flux}")
ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
ax2.set_ylabel(L"{\rm Flux}_{\rm GPU} - {\rm Flux}_{\rm CPU}")

ax1.legend()
fig.tight_layout()
plt.show()
