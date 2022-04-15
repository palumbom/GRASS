# imports
using Pkg; Pkg.activate(".")
using CUDA
using GRASS
using Statistics

# plotting imports
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# define some functions
include(GRASS.moddir * "figures1/fig_functions.jl")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# figure 2 -- cleaned + extrapolated input data
function main()
    # set number of curves to plot
    ncurves = 25

    # get list of directories
    dirs = filter(isdir, GRASS.soldir .* readdir(GRASS.soldir))

    for dir in dirs
        # figure out line name
        line_name = splitpath(dir)[end]
        line_title = replace(line_name, "_"=>" ")

        # get input data
        bisinfo = GRASS.SolarData(dir=dir, relative=true, extrapolate=true)
        key = (:c, :mu10)
        bis = bisinfo.bis[key]
        wav = bisinfo.wav[key]
        dep = bisinfo.dep[key]
        wid = bisinfo.wid[key]

        # create colormap and colorbar
        nstp = round(Int, size(bis,2) / ncurves)
        iter = 1:nstp:size(bis,2)
        iter = 1:2:30
        cmap = plt.get_cmap("Blues", length(iter))
        cols = [cmap(1*i/length(iter)) for i in 1:length(iter)]
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=15 * iter[end]/60.0))

        # plot the bisectors
        fig = plt.figure()
        ax1 = fig.add_subplot()
        for (i, t) in enumerate(iter)
            ax1.plot(wav[:,t][2:end].*1000, bis[:,t][2:end], c=cols[i], alpha=0.75)
        end

        # shade upper region
        xlims = ax1.get_xlim()
        ax1.fill_between(range(xlims..., step=0.1), 0.8, 1.0, hatch="/",
                         fc="black", ec="white", alpha=0.15, zorder=0)

        # set axis limits + labels
        ax1.set_xlim(xlims...)
        ax1.set_ylim(0.1, 1.0)
        ax1.set_title(latexstring(line_title))
        ax1.set_xlabel(L"{\rm Relative\ Wavelength\ (m\AA)}")
        ax1.set_ylabel(L"{\rm Normalized\ Intensity}")

        # set color bar
        cb = fig.colorbar(sm)
        cb.set_label(L"{\rm Time\ from\ first\ observation\ (min)}")

        # save the plot
        fig.savefig(plotdir * line_name * "_bisector.pdf")
        plt.clf(); plt.close()
        println(">>> Figure written to: " * plotdir * line_name * "_bisector.pdf")

        # plot the widths
        fig = plt.figure()
        ax1 = fig.add_subplot()
        for (i, t) in enumerate(iter)
            ax1.plot(dep[:,t], wid[:,t], c=cols[i], alpha=0.75)
        end

        # set axis limits + labels
        ax1.set_title(latexstring(line_title))
        ax1.set_xlabel(L"{\rm Normalized\ Intensity}")
        ax1.set_ylabel(L"{\rm Width\ across\ line\ (\AA)}")

        # set color bar
        cb = fig.colorbar(sm)
        cb.set_label(L"{\rm Time\ from\ first\ observation\ (min)}")

        # save the plot
        fig.savefig(plotdir * line_name * "_width.pdf")
        plt.clf(); plt.close()
        println(">>> Figure written to: " * plotdir * line_name * "_width.pdf")
    end
    return nothing
end

if (run | plot)
    main()
end
