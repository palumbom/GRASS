using Pkg; Pkg.activate(".")
using CSV
using GRASS
using Printf
using DataFrames
using Statistics
using EchelleCCFs

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use("my.mplstyle")

# set Desktop directory
outdir = "/Users/michael/Desktop/"
if !isdir(outdir)
    outdir = "/Users/mlp95/Desktop/"
end

# preprocess the data and get function definitions
function blueshift_vs_depth(depths::AbstractArray{Float64,1})
    # set up stuff for lines
    lines = [5434.5]
    resolution = 7e5
    templates = ["FeI_5434"]


    # allocate memory and loop over depths
    rvs = zeros(length(depths))
    for (ind, depth) in enumerate(depths)
        @printf("\t >>> Running depth %.2f \r", depth)
        spec = SpecParams(lines=lines, depths=[depth], templates=templates, resolution=resolution)

        # synthesize spectra
        disk = DiskParams(N=132, Nt=5)
        lambdas1, outspec1 = synthesize_spectra(spec, disk, verbose=false)

        # calculate the RVs
        v_grid, ccf1 = calc_ccf(lambdas1, outspec1, spec, normalize=true)
        rvs1, sigs1 = calc_rvs_from_ccf(v_grid, ccf1)
        rvs[ind] = mean(rvs1)
    end
    return rvs
end

# run it
depths = range(0.05, 0.95, step=0.1)
rvs = blueshift_vs_depth(depths)

# read in convective blueshifts from IAG paper
df_iag = GRASS.read_iag_blueshifts()

# get indices for bins
bin_edges = range(0.0, 1.0, step=0.1)
bin_centers = bin_edges[1:end-1] .+ 0.5 * step(bin_edges)
inds = zeros(Int, length(df_iag.depth))
for i in eachindex(df_iag.depth)
    inds[i] = findlast(x -> x .<= df_iag.depth[i], bin_edges)
end

# now bin the convective blueshifts
blueshift_avg = zeros(length(bin_centers))
blueshift_std = zeros(length(bin_centers))
blueshift_err = zeros(length(bin_centers))
for i in eachindex(bin_centers)
    data = filter(!isnan, df_iag.blueshift[inds .== i])
    blueshift_avg[i] = mean(data)
    blueshift_std[i] = std(data)
    blueshift_err[i] = blueshift_std[i] / sqrt(length(data))
end

# plot the result
fig, ax1 = plt.subplots()
ax1.scatter(df_iag.depth, df_iag.blueshift, s=1, alpha=0.5, c="black", label=L"{\rm IAG}", zorder=1)
ax1.errorbar(bin_centers, blueshift_avg, xerr=nothing, yerr=blueshift_err, fmt=".", capsize=3.0,
             markersize=10, alpha=0.85, c="black", label=L"{\rm Binned\ IAG}", zorder=1)
ax1.fill_between(bin_centers, blueshift_avg .- blueshift_std,
                 blueshift_avg .+ blueshift_std, color="k", alpha=0.25)
ax1.scatter(depths, rvs, c="tab:blue", s=15, alpha=0.85, label=L"{\rm Synthetic}", zorder=2)
ax1.set_ylim(-1100, 400)
ax1.set_xlabel(L"{\rm Line\ Depth}")
ax1.set_ylabel(L"{\rm Convective\ Blueshift\ (ms^{-1})}")
ax1.legend()
plt.show()
# fig.savefig(outdir * "conv_blueshift.pdf")
plt.clf(); plt.close()
