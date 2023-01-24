using Pkg; Pkg.activate(".")
using Statistics
using CUDA
using JLD2
using GRASS
using LsqFit
using FileIO
using DataFrames
using EchelleCCFs
using LaTeXStrings

# some global stuff
const N = round.(Int, 2 .^ range(6, 9, step=0.5))
const Nt = 100

# make data subdir
outdir = datadir * "resolutions/"
if !isdir(outdir)
    mkdir(outdir)
end

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

function resolution_curve(template; Nloop=16)
    # set up parameters for lines
    lines = [5434.5]
    depths = [0.8]
    templates = [template]
    res = 7e5
    spec = SpecParams(lines=lines, depths=depths, resolution=res, templates=templates)

    # allocate shared arrays
    avg_avg_res = zeros(length(N))
    std_avg_res = zeros(length(N))
    avg_rms_res = zeros(length(N))
    std_rms_res = zeros(length(N))

    # calculate
    @sync @distributed for i in 1:length(N)
    	println("running resolution N = " * string(N[i]))
        disk = DiskParams(N=N[i], Nt=Nt)
    	avg_avg1, std_avg1, avg_rms1, std_rms1 = GRASS.spec_loop(spec, disk, Nloop, use_gpu=true)
        avg_avg_res[i] = avg_avg1
        std_avg_res[i] = std_avg1
    	avg_rms_res[i] = avg_rms1
    	std_rms_res[i] = std_rms1
    end

    # make data frame
    df = DataFrame()
    df[!,:res] = N
    df[!,:avg_avg_res] = avg_avg_res
    df[!,:std_avg_res] = std_avg_res
    df[!,:avg_rms_res] = avg_rms_res
    df[!,:std_rms_res] = std_rms_res

    # write to CSV
    fname = datadir * "rms_vs_res.csv"
    CSV.write(fname, df)

    # write the results to file
    outfile = outdir * template * ".jld2"
    save(outfile,
         "res", res,
         "Nloop", Nloop,
         "template", template,
         "avg_avg_res", avg_avg_res,
         "std_avg_res", std_avg_res,
         "avg_rms_res", avg_rms_res,
         "std_rms_res", std_rms_res)
    return nothing
end

# run the simulation
if run
    # decide whether to use gpu
    @assert CUDA.functional()

    # loop over templates
    lp = GRASS.LineProperties()
    templates = GRASS.get_name(lp)
    for i in eachindex(templates)
        resolution_curve(templates[i])
    end
end

# plotting code block
if plot
    # plotting imports
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

    # make the plot subdir
    plotdir = plotdir * "resolutions/"
    if !isdir(plotdir)
        mkdir(plotdir)
    end

    # loop over line templates
    lp = GRASS.LineProperties()
    templates = GRASS.get_name(lp)
    for i in eachindex(templates)
        # get the filename and read in
        filename = outdir * templates[i] * ".jld2"

        d = load(filename)
        res = d["res"]
        Nloop = d["Nloop"]
        template = d["template"]
        rms_res = d["avg_rms_res"]
        std_res = d["std_rms_res"]

        # get the error
        err_res = std_res ./ sqrt(Nloop)

        # fit the data
        @. power_law(x, p) = p[1] * x^(-p[2])
        fit = curve_fit(power_law, res, rms_res, [1.0, 1.0])
        res_fit = range(6, 2400, length=1000)
        println(">>> Best fit power law index = " * string(fit.param[2]))
        x = string(round(fit.param[2], sigdigits=2))

        # initialize figure and set log scales
        fig = plt.figure()
        ax1 = fig.add_subplot()
        ax1.set_xscale("log", base=2)

        # plot the data
        ax1.errorbar(res, rms_res, yerr=err_res, capsize=3.0, color="black", fmt=".")
        ax1.plot(res_fit, power_law(res_fit, fit.param), "k--", alpha=0.4, label = L"{\rm Power\ law\ fit}")

        # add a line corresponding to average footprint area
        ax1.axvline(132, ymin=0, ymax=1, color="black", lw=1.5, ls="--")

        # plot the literature values
        xrng = ax1.get_xlim()
        xs = range(xrng[1], xrng[2], length=2)
        Elsworth = [repeat([0.319-0.09], 2), repeat([0.319+0.09], 2)]
        Palle = [repeat([0.461-0.1], 2), repeat([0.461+0.1], 2)]
        ax1.fill_between(xs, Elsworth..., alpha=0.5, fc="tab:orange", ec="white",
                         hatch="/", label=L"{\rm Elsworth\ et\ al.\ (1994)}")
        ax1.fill_between(xs, Palle..., alpha=0.5, fc="tab:green", ec="white",
                         hatch="\\", label=L"\textnormal{Pall\'{e}}\ {\rm et\ al.\ (1999)}")

        # set the axes limits + labels
        ax1.set_xlim(2^5.9, 2^10.1)
        ax1.set_ylim(0.01, 1.15)
        ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

        # set titles
        txt = ("\${\\rm " * replace(template, "_" => "\\ ") * "}\$")
        ax1.set_title(txt)
        ax1.set_xlabel(L"N")
        ax1.set_ylabel(L"{\rm RMS}_{\rm RV}\ {\rm (m s}^{-1})")
        ax1.legend()

        # save the figure
        fig.savefig(plotdir * template * "_res.pdf")
        plt.clf(); plt.close()
        println(">>> Figure written to: " * plotdir * template * "_res.pdf")
    end
end
