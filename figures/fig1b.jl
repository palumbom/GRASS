# imports
using Pkg; Pkg.activate(".")
using GRASS
using Statistics

# plotting imports
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

# define some functions
include(GRASS.moddir * "figures/fig_functions.jl")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# figure 1b -- input bisectors w/ variability
function main()
    # get input data
    bisinfo = GRASS.SolarData()

    # initialize plot objects
    fig, ax1 = plt.subplots()

    # loop and plot curves
    keyz = [(:c, :mu10), (:w, :mu06), (:w, :mu03)]
    labels = [L"\mu = 1.0", L"\mu = 0.6", L"\mu = 0.3"]
    linestyles = ["-", "--", ":"]
    for (i, key) in enumerate(keyz)
        # find average and std
        avg_bis = mean(bisinfo.bis[key], dims=2)
        avg_wav = mean(bisinfo.wav[key], dims=2)
        std_bis = std(bisinfo.bis[key], dims=2)
        std_wav = std(bisinfo.wav[key], dims=2)

        # convert to doppler velocity
        restwav = 5434.5232                         # angstroms
        avg_wav = avg_wav ./ restwav .* GRASS.c_ms
        std_wav = std_wav ./ restwav .* GRASS.c_ms

        # cut off top portion, where uncertainty is large
        ind = findfirst(avg_bis .> 0.86)[1]
        avg_bis = avg_bis[1:ind]
        avg_wav = avg_wav[1:ind]
        std_bis = std_bis[1:ind]
        std_wav = std_wav[1:ind]

        # fix dimensions
        y = reshape(avg_bis, length(avg_bis))
        x1 = reshape(avg_wav .+ std_wav, length(avg_bis))
        x2 = reshape(avg_wav .- std_wav, length(avg_bis))

        # plot the curve
        ax1.fill_betweenx(y, x1, x2, color="C"*string(i-1), alpha=0.5)
        ax1.plot(avg_wav, avg_bis, ls=linestyles[i], color="C"*string(i-1), label=labels[i])
    end

    # set axes labels and save the figure
    ax1.legend(loc="upper right")
    ax1.set_ylim()
    ax1.set_xlabel(L"{\rm Doppler\ Velocity\ (ms}^{-1} {\rm )}")
    ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
    fig.savefig(plotdir * "fig1b.pdf")
    plt.clf(); plt.close()
    println(">>> Figure written to: " * plotdir * "fig1b.pdf")
    return nothing
end

if (run | plot)
    main()
end
