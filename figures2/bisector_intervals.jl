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

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# create output directory if needed
outdir = joinpath(plotdir, "bisector_intervals")
if !isdir(outdir)
    mkdir(outdir)
end

# figure 1b -- input bisectors w/ variability
function main()
    # looping over inputs
    lp = LineProperties()
    files = lp.file
    names = GRASS.get_name(lp)
    restwavs = GRASS.get_rest_wavelength(lp)

    for (idx, fname) in enumerate(files)
        # figure out line name
        line_name = names[idx]
        line_title = replace(line_name, "_"=>" ")

        # get input data
        bisinfo = GRASS.SolarData(fname=fname)

        # initialize plot objects

        fig, ax1 = plt.subplots()

        # loop and plot curves
        keyz = [(:c, :mu10), (:w, :mu06), (:w, :mu03)]
        labels = [L"\mu = 1.0", L"\mu = 0.6", L"\mu = 0.3"]
        linestyles = ["-", "--", ":"]
        for (i, key) in enumerate(keyz)
            if !(key in keys(bisinfo.bis))
                continue
            end

            # find average and std
            avg_bis = dropdims(mean(bisinfo.bis[key], dims=2), dims=2)
            avg_int = dropdims(mean(bisinfo.int[key], dims=2), dims=2)
            std_bis = dropdims(std(bisinfo.bis[key], dims=2), dims=2)
            std_int = dropdims(std(bisinfo.int[key], dims=2), dims=2)

            # convert to doppler velocity
            avg_bis .*= GRASS.c_ms / (restwavs[idx])
            std_bis .*= GRASS.c_ms / (restwavs[idx])

            # cut off top portion, where uncertainty is large
            # ind = findfirst(avg_int .> 0.86)[1]
            # avg_bis = avg_bis[1:ind]
            # avg_int = avg_int[1:ind]
            # std_bis = std_bis[1:ind]
            # std_int = std_int[1:ind]

            # fix dimensions
            y = reshape(avg_int, length(avg_int))
            x1 = reshape(avg_bis .+ std_bis, length(avg_int))
            x2 = reshape(avg_bis .- std_bis, length(avg_int))

            # plot the curve
            ax1.fill_betweenx(y, x1, x2, color="C"*string(i-1), alpha=0.5)
            ax1.plot(avg_bis, avg_int, ls=linestyles[i], color="C"*string(i-1), label=labels[i])
        end

        # set axes labels and save the figure
        ax1.legend(loc="upper right")
        ax1.set_xlim(-250, 250)
        # ax1.set_ylim()
        ax1.set_title("\${\\rm " * replace(line_name, "_" => "\\ ") * "}\$")
        ax1.set_xlabel(L"{\rm Doppler\ Velocity\ (ms}^{-1} {\rm )}")
        ax1.set_ylabel(L"{\rm Normalized\ Intensity}")

        # save it
        outname = joinpath(outdir, line_name * "_interval.pdf")
        fig.savefig(outname)
        plt.clf(); plt.close()
        println(">>> Figure written to: " * outname)
    end
    return nothing
end

if (run | plot)
    main()
end
