# imports
using Pkg; Pkg.activate(".")
using GRASS
using Statistics

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

# define rms loop function
include(GRASS.moddir * "figures/rms_loop.jl")

# set boolean for writing plot
write = true
grassdir, plotdir, datadir = check_plot_dirs()

# figure 2 -- cleaned + extrapolated input data
function plot_input_cleaned(ncurves)
    # get input data
    bisinfo = GRASS.SolarData(relative=true, extrapolate=true)
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
        ax1.plot(wav[:,t].*1000, bis[:,t], c=cols[i], alpha=0.75)
    end
    ax1.set_xlabel(L"{\rm Relative\ Wavelength\ (m\AA)}")
    ax1.set_ylabel(L"{\rm Normalized\ Flux}")
    cb = fig.colorbar(sm)
    cb.set_label(L"{\rm Time\ from\ first\ observation\ (min)}")

    # write the file or show it
    if write
        fig.savefig(plotdir * "fig2a.pdf"))
        plt.clf(); plt.close()
    else
        plt.show()
        plt.clf(); plt.close()
    end

    # plot the widths
    fig = plt.figure()
    ax1 = fig.add_subplot()
    for (i, t) in enumerate(iter)
        ax1.plot(dep[:,t], wid[:,t], c=cols[i], alpha=0.75)
    end
    ax1.set_xlabel(L"{\rm Normalized\ Flux}")
    ax1.set_ylabel(L"{\rm Width\ across\ line\ (\AA)}")
    cb = fig.colorbar(sm)
    cb.set_label(L"{\rm Time\ from\ first\ observation\ (min)}")

    # write the file or show it
    if write
        fig.savefig(plotdir * "fig2b.pdf"))
        plt.clf(); plt.close()
    else
        plt.show()
        plt.clf(); plt.close()
    end
    return nothing
end

plot_input_cleaned(25)
