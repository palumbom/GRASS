# environment + packages
using Pkg; Pkg.activate(".")
using CSV
using HTTP
using GZip
using GRASS
using LsqFit
using Statistics
using DataFrames
using Interpolations

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

# # define some functions
include(GRASS.moddir * "figures/fig_functions.jl")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

function download_iag()
    println(">>> Downloading IAG atlas...")
    file = HTTP.download("https://cdsarc.unistra.fr/ftp/J/A+A/587/A65/spvis.dat.gz",
                         GRASS.moddir * "figures/", update_period=Inf)
    println(">>> IAG atlas downloaded!")
    return nothing
end

function read_iag(; isolate=true)
    # download the IAG atlas
    file = GRASS.moddir * "figures/spvis.dat.gz"
    if !isfile(file)
        download_iag()
    end

    # read in the IAG atlas
    iag = GZip.open(file, "r") do io
        CSV.read(io, DataFrame, ignorerepeated=true, delim=" ", header=["wavenum", "nflux", "flux"])
    end

    # convert wavenumber to wavelength in angstroms
    wavs = (1 ./ iag.wavenum) * 1e8

    # reverse to deal with conversion of units
    reverse!(wavs)
    reverse!(iag.nflux)

    # isolate region around 5434.5 line
    if isolate
        ind1 = findfirst(x -> x .> 5435.5, wavs)
        ind2 = findfirst(x -> x .> 5436.55, wavs)
        return wavs[ind1:ind2], iag.nflux[ind1:ind2]
    end
    return wavs, iag.nflux[ind1:ind2]
end

function interpolate_spec(wavs, flux)
    wavs_itp = collect(range(wavs[1], wavs[end], step=mean(diff(wavs))))
    itp = LinearInterpolation(wavs, flux)
    flux_itp = itp.(wavs_itp)
    return wavs_itp, flux_itp
end

# model the iag blends
function model_iag_blends(wavs_sim, flux_sim, wavs_iag, flux_iag; plot=false)
    # roughly align the spectra
    off1 = wavs_iag[argmin(flux_iag)] - wavs_sim[argmin(flux_sim)]
    wavs_iag .-= off1

    # interpolate onto the same grid
    flux_iag ./= maximum(flux_iag)
    itp = LinearInterpolation(wavs_sim, flux_sim, extrapolation_bc=1.0)
    wavs_sim = wavs_iag
    flux_sim = itp.(wavs_sim)

    # calculate some residuals
    resids = flux_iag .- flux_sim

    # models for fit
    @. gaussian(x, a, b, c) = a * exp(-(x - b)^2/(2 * c^2)) + 1
    function tel_model(x, p)
        n = length(p) ÷ 3
        out = ones(length(x))
        for i in 1:n
            out .*= gaussian(x, p[3i-2:3i]...)
        end
        return out .* flux_sim
    end
    # initial guesses for blended lines
    l0 = [-0.005, 5434.0, 0.012]
    l1 = [-0.03, 5434.15, 0.0325]
    l2 = [-0.02, 5434.3, 0.03]
    l3 = [-0.03, 5434.48, 0.02]
    l4 = [-0.02, 5434.55, 0.025]
    l5 = [-0.03, 5434.65, 0.06]
    l6 = [-0.025, 5434.71, 0.02]
    l7 = [-0.05, 5434.812, 0.0305]
    l8 = [-0.035, 5435.01, 0.04]
    l9 = [-0.01, 5435.025, 0.02]
    p0 = [l1..., l2..., l4..., l6..., l7..., l8...]

    # do the fit
    fit = curve_fit(tel_model, wavs_iag, flux_iag, p0)

    # plot diagnostics
    if plot
        fig = plt.figure(figsize=(8,6))
        gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
        ax1 = fig.add_subplot(gs[1])
        ax2 = fig.add_subplot(gs[2])
        ax1.plot(wavs_iag, flux_iag, label="Observed IAG", color="tab:blue")
        ax1.plot(wavs_iag, tel_model(wavs_iag, fit.param), label="Modeled IAG ", ls=":", color="tab:green")
        ax2.scatter(wavs_sim, flux_iag./tel_model(wavs_iag, fit.param), c="k", s=0.5)
        ax1.legend()
        ax1.set_xticklabels([])
        ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
        ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
        ax2.set_ylabel(L"{\rm IAG/Model}")
        plt.show()
        plt.clf(); plt.close()
    end
    return flux_iag./(tel_model(wavs_iag, fit.param)./flux_sim)
end

# figure 3 -- compare synthetic and IAG spectra + bisectors
function main()
    # some parameters
    vlim = 2.1e4
    btop = 0.9

    # get spectrum and interpolate onto even wavelength grid
    wavs_iag, flux_iag = read_iag(isolate=true)
    wavs_iag, flux_iag = interpolate_spec(wavs_iag, flux_iag)

    # set up for GRASS spectrum simulation
    lines = [5434.5232]
    depths = [1.0 - minimum(flux_iag)]
    resolution = 700000.0
    spec = SpecParams(lines=lines, depths=depths, resolution=resolution, extrapolate=true)
    disk = DiskParams(N=132, Nt=15)

    # synthesize spectra, calculate ccf, and get CCF bisector
    len = 72
    lambdas1, outspec1 = synthesize_spectra(spec, disk, seed_rng=false, verbose=true, top=NaN)
    outspec1 ./= maximum(outspec1)
    v_grid, ccf1 = calc_ccf(lambdas1, outspec1, spec, normalize=true)
    outspec1 = mean(outspec1, dims=2)[:,1]
    ccfm = mean(ccf1, dims=2)[:,1]
    vel_sim, bis_sim = GRASS.measure_bisector(v_grid, ccfm, interpolate=false, top=btop, len=len)

    # model the line blends out of the IAG spectrum
    println(">>> Modeling out line blends in IAG spectrum...")
    flux_iag_cor = model_iag_blends(lambdas1, outspec1, wavs_iag, flux_iag, plot=false)

    # get offsets to align the spectra in wavelength
    off1 = wavs_iag[argmin(flux_iag)] - lambdas1[argmin(outspec1[:,1]), 1]
    off2 = wavs_iag[argmin(flux_iag_cor)+1] - lambdas1[argmin(outspec1[:,1]), 1]
    wavs_iag .-= off1

    # calculate a CCF for the IAG spectrum and trim it
    v_grid_iag, ccf_iag = calc_ccf(wavs_iag, flux_iag, [wavs_iag[argmin(flux_iag)]],
                                   [1.0 - minimum(flux_iag)],
                                   1e6, normalize=true)
    ind1 = findfirst(x -> x .> -vlim, v_grid_iag)
    ind2 = findfirst(x -> x .> vlim, v_grid_iag)
    vel_iag, bis_iag = GRASS.measure_bisector(v_grid_iag[ind1:ind2], ccf_iag[ind1:ind2],
                                              interpolate=true, top=btop, len=len)

    # calculate a CCF for the cleaned IAG spectrum and trim it
    wavs_iag, flux_iag_cor = interpolate_spec(wavs_iag, flux_iag_cor)
    v_grid_iag2, ccf_iag2 = calc_ccf(wavs_iag, flux_iag_cor, [wavs_iag[argmin(flux_iag_cor)]],
                                    [1.0 - minimum(flux_iag_cor)],
                                    1e6, normalize=true)
    ind1 = findfirst(x -> x .> -vlim, v_grid_iag2)
    ind2 = findfirst(x -> x .> vlim, v_grid_iag2)
    vel_iag2, bis_iag2 = GRASS.measure_bisector(v_grid_iag2[ind1:ind2], ccf_iag2[ind1:ind2],
                                                interpolate=true, top=btop, len=len)

    # interpolate IAG onto same wavelength scale as synthetic spectrum
    itp = LinearInterpolation(wavs_iag, flux_iag, extrapolation_bc=1.0)
    flux_iag_itp = itp.(lambdas1)
    itp = LinearInterpolation(wavs_iag, flux_iag_cor, extrapolation_bc=1.0)
    flux_iag_itp2 = itp.(lambdas1)

    function comparison_plots()
        # overplot the spectra
        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
        ax1 = fig.add_subplot(gs[1])
        ax2 = fig.add_subplot(gs[2])
        ax1.plot(lambdas1, outspec1, c="black", lw= 1.5, label=L"{\rm Synthetic}")
        ax1.plot(wavs_iag, flux_iag./maximum(flux_iag), marker="s", c="tab:blue", ms=2.0, lw=1.0, markevery=10, label=L"{\rm IAG}")
        ax1.plot(wavs_iag, flux_iag_cor./maximum(flux_iag_cor), alpha=0.9, marker="o", c="tab:green", ms=2.0, lw=1.0, markevery=10, label=L"{\rm Cleaned\ IAG}")
        ax2.plot(lambdas1, flux_iag_itp./maximum(flux_iag_itp) .- outspec1, c="tab:blue", marker="s", ms=2.0, lw=0)
        ax2.plot(lambdas1, flux_iag_itp2./maximum(flux_iag_itp2) .- outspec1, c="tab:green", marker="o", ms=2.0, lw=0)

        # set tick labels, axis labels, etc.
        ax1.set_xticklabels([])
        ax1.set_xlim(5434, 5435)
        ax2.set_xlim(5434, 5435)
        ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
        ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
        ax2.set_ylabel(L"{\rm IAG\ -\ Synthetic}")
        ax1.legend()
        fig.tight_layout()

        # save the plot
        fig.savefig(plotdir * "fig3a.pdf")
        plt.clf(); plt.close()
        println(">>> Figure written to: " * plotdir * "fig3a.pdf")

        # align bisectors to arbitrary point
        vel_sim .-= mean(vel_sim)
        vel_iag .-= mean(vel_iag) - 42.0
        vel_iag2 .-= mean(vel_iag2)

        # plot the bisectors
        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
        ax1 = fig.add_subplot(gs[1])
        ax2 = fig.add_subplot(gs[2])
        ax1.plot(vel_sim, bis_sim, color="black", lw=2.0, label=L"{\rm Synthetic}")
        ax1.plot(vel_iag, bis_iag, marker="s", c="tab:blue", ms=2.0, lw=1.0, label=L"{\rm IAG}")
        ax1.plot(vel_iag2, bis_iag2,marker="o", c="tab:green", ms=2.0, lw=1.0, label=L"{\rm Cleaned\ IAG}")
        ax2.plot(vel_iag .- vel_sim, bis_iag, c="tab:blue", marker="s", ms=2.0, lw=0.0)
        ax2.plot(vel_iag2 .- vel_sim, bis_iag, c="tab:green", marker="o", ms=2.0, lw=0.0)

        # set tick labels, axis labels, etc.
        ax2.set_yticklabels([])
        ax2.yaxis.tick_right()
        ax1.set_ylim(0.1, 1.1)
        ax2.set_xlim(-15, 25)
        ax2.set_ylim(0.1, 1.1)
        ax1.set_xlabel(L"{\rm Relative\ Velocity\ (ms^{-1})}")
        ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
        ax2.set_xlabel(L"{\rm IAG\ -\ Synthetic\ (ms^{-1})}")
        ax1.legend(loc="upper right", prop=Dict("size"=>10), labelspacing=0.25)

        # save the plot
        fig.savefig(plotdir * "fig3b.pdf")
        plt.clf(); plt.close()
        println(">>> Figure written to: " * plotdir * "fig3b.pdf")
        return nothing
    end

    comparison_plots()
    return nothing
end

if (run | plot)
    main()
end
