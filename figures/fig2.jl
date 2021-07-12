using Pkg; Pkg.activate(".")
using Profile
using Statistics
using BenchmarkTools
using Printf

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use("my.mplstyle")

# SynthSpectra
using Revise
using SynthSpectra; SS=SynthSpectra
using EchelleCCFs
using Debugger

# set Desktop directory
outdir = "/Users/michael/Desktop/"
if !isdir(outdir)
    outdir = "/Users/mlp95/Desktop/"
end

# figure 1 -- input spectra w/ variability
function plot_input_variability()
    # get input data
    bisinfo = SS.SolarData(relative=true)

    # loop and plot
    keyz = [(:c, :mu10), (:w, :mu06), (:w, :mu03)]
    labels = [L"\mu = 1.0", L"\mu = 0.6", L"\mu = 0.3"]
    for (i, key) in enumerate(keyz)
        # find average and std
        avg_bis = mean(bisinfo.bis[key], dims=2)
        avg_wav = mean(bisinfo.wav[key], dims=2)
        std_bis = std(bisinfo.bis[key], dims=2)
        std_wav = std(bisinfo.wav[key], dims=2)

        # convert to doppler velocity
        avg_wav = avg_wav ./ 5434.5232 .* SS.c_ms
        std_wav = std_wav ./ 5434.5232 .* SS.c_ms

        # cut off top
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
        plt.fill_betweenx(y, x1, x2, color="C"*string(i-1), alpha=0.5)
        plt.plot(avg_wav, avg_bis, color="C"*string(i-1), label=labels[i])
    end
    plt.legend(loc="upper right")
    plt.xlabel(L"{\rm Doppler\ Velocity\ (ms}^{-1} {\rm )}")
    plt.ylabel(L"{\rm Normalized\ Flux}")
    plt.savefig(outdir * "fig1b.pdf", bbox_inches="tight")
    plt.close(); plt.clf()
    return nothing
end

plot_input_variability()

# figure 2 -- input data
function plot_input_cleaned(ncurves)
    # get input data
    bisinfo = SS.SolarData(relative=true, extrapolate=true)
    key = (:c, :mu10)
    bis = bisinfo.bis[key]
    wav = bisinfo.wav[key]
    dep = bisinfo.dep[key]
    wid = bisinfo.wid[key]

    # cmap stuff
    nstp = round(Int, size(bis,2) / ncurves)
    iter = 1:nstp:size(bis,2)
    iter = 1:2:30
    cmap = plt.get_cmap("Blues", length(iter))
    cols = [cmap(1*i/length(iter)) for i in 1:length(iter)]

    # colorbar stuff
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
    fig.savefig(outdir * "fig2a.pdf")
    plt.clf(); plt.close()

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
    fig.savefig(outdir * "fig2b.pdf")
    plt.clf(); plt.close()
    return nothing
end

plot_input_cleaned(25)


# bisector span
# function plot_bisector_span()
#     bisinfo = SS.BisectorInfo(relative=false)

#     # loop and plot
#     keyz = keys(bisinfo.wav)
#     avg_span = []
#     std_span = []
#     mu = []
#     mu_symb = [:mu02, :mu03, :mu04, :mu05, :mu06, :mu07, :mu08, :mu085, :mu09, :mu095, :mu10]
#     disc_mu = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0]
#     for (i, key) in enumerate(keys(bisinfo.wav))
#         !(:e in key) && continue
#         mukey = key[2]
#         ind = findfirst(mukey .== mu_symb)

#         spans = SS.calculate_bisector_span(bisinfo.line.λrest, bisinfo.wav[key])
#         push!(avg_span, mean(spans))
#         push!(std_span, std(spans))
#         push!(mu, disc_mu[ind])
#     end
#     spans = SS.calculate_bisector_span(bisinfo.line.λrest, bisinfo.wav[(:c, :mu10)])
#     push!(avg_span, mean(spans))
#     push!(std_span, std(spans))
#     push!(mu, 1.0)

#     fig = plt.figure()
#     ax1 = fig.add_subplot()
#     ax1.errorbar(mu, avg_span, yerr=std_span, capsize=3.0, fmt="k.")
#     ax1.ticklabel_format(useOffset=false)
#     plt.locator_params(nbins=5,axis="x")
#     ax1.set_xlabel(L"\mu")
#     ax1.set_ylabel(L"{\rm Bisector\ Span}\ ({\rm ms}^{-1})")
#     ax1.invert_xaxis()
#     fig.savefig(outdir * "fig2.pdf")
#     plt.clf(); plt.close()
# end

