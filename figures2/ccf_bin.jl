using Pkg; Pkg.activate(".")
using CSV
using CUDA
using JLD2
using GRASS
using Printf
using Revise
using FileIO
using DataFrames
using Statistics
using EchelleCCFs
using Polynomials
using Distributions
using BenchmarkTools
using HypothesisTests

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")

# # get input data properties
# lp = GRASS.LineProperties()
# line_species = GRASS.get_species(lp)
# rest_wavelengths = GRASS.get_rest_wavelength(lp)
# line_depths = GRASS.get_depth(lp)
# line_names = GRASS.get_name(lp)
# line_titles = replace.(line_names, "_" => " ")
# line_files = GRASS.get_file(lp)

# # draw random line depths and centers
# nlines = 10
# lines = Float64[]
# depths = Float64[]
# templates = String[]
# for i in eachindex(rest_wavelengths)
#     # ltemp = rand(Uniform(minimum(rest_wavelengths), maximum(rest_wavelengths)), nlines)
#     ltemp = rand(Uniform(5200, 5400), nlines)
#     dtemp = rand(Normal(line_depths[i]-0.05, 0.05), nlines)
#     push!(lines, ltemp...)
#     push!(depths, dtemp...)
#     push!(templates, repeat([line_files[i]], nlines)...)
# end

# # synthesize a spectrum
# N = 132
# Nt = 1000
# variability = trues(length(lines))
# resolution = 7e5
# seed_rng = true

# disk = DiskParams(N=N, Nt=Nt)
# spec1 = SpecParams(lines=lines, depths=depths, variability=variability, templates=templates, resolution=resolution)
# lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=true, verbose=true, use_gpu=true)

# # save it to a JLD
# jldsave("spectra_for_bin.jld2", wavs=lambdas1, flux=outspec1, templates=templates, lines=lines, depths=depths)

function make_ccf_plots(wavs, flux, lines, depths, title)
    # calculate a ccf for one line
    v_grid, ccf1 = calc_ccf(wavs, flux, lines, depths,
                            7e5, mask_type=EchelleCCFs.TopHatCCFMask,
                            Δv_step=125.0)
    rvs1, sigs1 = calc_rvs_from_ccf(v_grid, ccf1)
    bis, int = GRASS.calc_bisector(v_grid, ccf1, nflux=20)
    bis_inv_slope = GRASS.calc_bisector_inverse_slope(bis, int)

    # get data to fit / plot
    xdata = rvs1 .- mean(rvs1)
    ydata = bis_inv_slope .- mean(bis_inv_slope)

    # do the fit
    pfit = Polynomials.fit(xdata, ydata, 1)
    xmodel = range(-2, 2, length=5)
    ymodel = pfit.(xmodel)

    # get slope string to annotate
    slope = round(coeffs(pfit)[2], digits=3)
    fit_label = "\$ " .* string(slope) .* "\$"

    # calc mean to plot
    mean_bis = dropdims(mean(bis, dims=2), dims=2)[2:end-2]
    mean_int = dropdims(mean(int, dims=2), dims=2)[2:end-2]
    mean_bis .-= mean(mean_bis)

    # make plotting objects
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 7))
    fig.subplots_adjust(wspace=0.05)
    # ax1.set_box_aspect(1)
    # ax2.set_box_aspect(1)

    ax1.plot(mean_bis, mean_int, c="k", lw=3.0)
    ax2.scatter(xdata, ydata, c="k", s=1.5, alpha=0.9)
    ax2.plot(xmodel, ymodel, c="k", ls="--", label=L"{\rm Slope } \approx\ " * fit_label, lw=2.5)
    ax2.legend()

    ax1.set_xlabel(L"\Delta\ v\ {\rm (m/s)}")
    ax1.set_ylabel(L"{\rm Normalized\ CCF}")
    ax2.set_xlabel(L"\Delta v\ {\rm (m s}^{-1}{\rm )}")
    ax2.set_ylabel(L"{\rm BIS}\ - \overline{\rm BIS}\ {\rm (m s}^{-1}{\rm )}")
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    fig.suptitle("{\\rm " * replace(title, " "=> "\\ ") * "}")
    fig.tight_layout()
    fig.savefig(joinpath("/Users/michael/Desktop/ccf_bin_plots/", title * ".pdf"))
    # plt.show()
    plt.clf(); plt.close()
end

# load the file
data = load("/Users/michael/Desktop/spectra_for_bin.jld2")
depths = data["depths"]
wavs = data["wavs"]
lines = data["lines"]
flux = data["flux"]
templates = data["templates"]

# do a random line and do all lines
make_ccf_plots(wavs, flux, [lines[25]], [depths[25]], "Single Line")
make_ccf_plots(wavs, flux, lines, depths, "All Lines")

# get list of line idxs in each depth bin
depth_bins = range(0.0, 1.0, step=0.1)
idxs = []
for i in eachindex(depth_bins)
    # get indices
    i == 1 && continue
    idx = map(x -> (x .<= depth_bins[i]) & (x .> depth_bins[i - 1]), depths)
    push!(idxs, idx)
end

#
for i in eachindex(idxs)
    # get line title
    title = "Depths gt " * string(depth_bins[i]) * " and lt " * string(depth_bins[i+1])
    println(title)

    # get the lines and make plots
    lines_i = view(lines, idxs[i])
    depths_i = view(depths, idxs[i])
    if isempty(lines_i)
        continue
    end
    make_ccf_plots(wavs, flux, lines_i, depths_i, title)
end

# read in formation temps
line_info = CSV.read(GRASS.datdir * "line_info.csv", DataFrame)
