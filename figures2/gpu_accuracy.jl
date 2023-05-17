using Pkg; Pkg.activate(".")
using JLD2
using CUDA
using GRASS
using Printf
using FileIO
using Revise
using Formatting
using Statistics
using EchelleCCFs
using Distributions
using BenchmarkTools
using HypothesisTests

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

# parse args + get directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

"""
# set up paramaters for spectrum
lines = [5500.0, 5500.85, 5501.4, 5502.20, 5502.5, 5503.05]
depths = [0.75, 0.4, 0.65, 0.55, 0.25, 0.7]
templates = ["FeI_5434", "FeI_5576", "FeI_6173", "FeI_5432", "FeI_5576", "FeI_5250.6"]
resolution = 7e5
buffer = 0.6

# create composite types
disk = DiskParams(N=132, Nt=5)
spec = SpecParams(lines=lines, depths=depths, templates=templates,
                  resolution=resolution, buffer=buffer)

# get the spectra
println(">>> Doing CPU synthesis...")
wavs_cpu64, flux_cpu64 = synthesize_spectra(spec, disk, seed_rng=true, verbose=true, use_gpu=false)
println(">>> Doing GPU synthesis (double precision)...")
wavs_gpu64, flux_gpu64 = synthesize_spectra(spec, disk, seed_rng=true, verbose=true, use_gpu=true)
println(">>> Doing GPU synthesis (single precision)...")
wavs_gpu32, flux_gpu32 = synthesize_spectra(spec, disk, seed_rng=true, verbose=true, use_gpu=true, precision=Float32)

# compute means
flux_cpu_mean64 = dropdims(mean(flux_cpu64, dims=2), dims=2)
flux_gpu_mean64 = dropdims(mean(flux_gpu64, dims=2), dims=2)
flux_gpu_mean32 = dropdims(mean(flux_gpu32, dims=2), dims=2)

# get residuals
resids64 = flux_cpu_mean64 .- flux_gpu_mean64
resids32 = flux_cpu_mean64 .- flux_gpu_mean32

# test residuals for normality
AD1 = OneSampleADTest(resids64, Normal())
AD2 = OneSampleADTest(resids32, Normal())
"""

# set up plot
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(8,10), sharex=true)

ms = 2.0
ax1.plot(wavs_cpu64, flux_cpu_mean64, ls="--", c="k", ms=ms, label=L"{\rm CPU\ (Float64)}")
ax1.plot(wavs_gpu64, flux_gpu_mean64, ls="-.", c=colors[1], ms=ms, label=L"{\rm GPU\ (Float64)}")
ax1.plot(wavs_gpu32, flux_gpu_mean32, ls=":", c=colors[2], ms=ms, label=L"{\rm GPU\ (Float32)}")
ax2.scatter(wavs_cpu64, resids64, s=5, marker="o", c=colors[1], alpha=0.9)
ax3.scatter(wavs_cpu64, resids32, s=5, marker="s", c=colors[2], alpha=0.9)

# inset axis for marginal distribution of residuals
ax2_histy = ax2.inset_axes([1.0015, 0, 0.075, 1], sharey=ax2)
ax2_histy.hist(resids64, bins="auto", density=true, histtype="step", orientation="horizontal", color=colors[1])
ax2_histy.tick_params(axis="both", labelleft=false, labelbottom=false)
ax2_histy.set_xticks([])
ax2_histy.get_yaxis().set_visible(false)
ax2_histy.spines["left"].set_visible(false)
ax2_histy.spines["bottom"].set_visible(false)
ax2_histy.spines["top"].set_visible(false)
ax2_histy.spines["right"].set_visible(false)
ax2_histy.grid(false)

ax3_histy = ax3.inset_axes([1.0015, 0, 0.075, 1], sharey=ax3)
ax3_histy.hist(resids32, bins="auto", density=true, histtype="step", orientation="horizontal", color=colors[2])
ax3_histy.tick_params(axis="both", labelleft=false, labelbottom=false)
ax3_histy.set_xticks([])
ax3_histy.get_yaxis().set_visible(false)
ax3_histy.spines["left"].set_visible(false)
ax3_histy.spines["bottom"].set_visible(false)
ax3_histy.spines["top"].set_visible(false)
ax3_histy.spines["right"].set_visible(false)
ax3_histy.grid(false)

ax1.set_xlim(minimum(wavs_cpu64)-0.1, maximum(wavs_cpu64)+0.1)
ax1.set_ylim(0.075, 1.05)
ax2.set_ylim(-1.2e-14, 1.2e-14)
ax3.set_ylim(-1.2e-3, 1.2e-3)

# do manual formatting of tick labels since the default is ROUGH
ax2.ticklabel_format(axis="y", useOffset=true, style="sci")
ax2.yaxis.offsetText.set_visible(false)

ax3.ticklabel_format(axis="y", useOffset=true, style="sci")
ax3.yaxis.offsetText.set_visible(false)
ax3.set_yticks([-1.0e-3, -0.5e-3, 0.0, 0.5e-3, 1.0e-3])
ax3.set_yticklabels([L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"])

ax1.legend(loc="lower center", ncols=3, fontsize=13)
ax1.set_ylabel(L"{\rm Normalized\ Flux}")
ax2.set_ylabel(L"({\rm Flux}_{\rm CPU} - {\rm Flux}_{\rm GPU}) \times 10^{-15}")
ax3.set_ylabel(L"({\rm Flux}_{\rm CPU} - {\rm Flux}_{\rm GPU}) \times 10^{-3}")
ax3.set_xlabel(L"{\rm Wavelength\ (\AA)}")

fig.tight_layout()
fig.subplots_adjust(hspace=0.05)
fig.savefig(joinpath(plotdir, "gpu_accuracy.pdf"))
plt.clf(); plt.close()
