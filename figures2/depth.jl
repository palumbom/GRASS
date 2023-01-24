# import stuff
using Pkg; Pkg.activate(".")
using CUDA
using JLD2
using GRASS
using FileIO
using DataFrames
using LaTeXStrings

# some global stuff
const N = 132
const Nt = 300

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# make data subdir
outdir = datadir * "depths/"
if !isdir(outdir)
    mkdir(outdir)
end

function single_line_variability(wavelength::Float64, template::String; Nloop::Int=25)
    # set up parameters for lines
    lines = [wavelength]
    depths = range(0.05, stop=0.95, step=0.05)
    resolution=7e5
    use_gpu = true

    # allocate shared arrays
    avg_avg_depth = zeros(length(depths))
    std_avg_depth = zeros(length(depths))
    avg_rms_depth = zeros(length(depths))
    std_rms_depth = zeros(length(depths))

    # loop over depths
    disk = DiskParams(N=N, Nt=Nt)
    for i in eachindex(depths)
        println("running depth = " * string(depths[i]))

        # create spec instance
        spec = SpecParams(lines=lines, depths=[depths[i]], templates=[template], resolution=resolution)

        # synthesize spectra, get velocities and stats
        out = GRASS.spec_loop(spec, disk, Nloop, use_gpu=use_gpu)
        avg_avg_depth[i] = out[1]
        std_avg_depth[i] = out[2]
        avg_rms_depth[i] = out[3]
        std_rms_depth[i] = out[4]
    end

    # write the results to file
    outfile = outdir * template * ".jld2"
    save(outfile,
         "Nloop", Nloop,
         "wavelength", wavelength,
         "template", template,
         "depths", depths,
         "avg_avg_depth", avg_avg_depth,
         "std_avg_depth", std_avg_depth,
         "avg_rms_depth", avg_rms_depth,
         "std_rms_depth", std_rms_depth)
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
        single_line_variability(lp.λrest[i], templates[i])
    end
end

if plot
    # plotting imports
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

    # make the plot subdir
    plotdir = plotdir * "depths/"
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
        Nloop = d["Nloop"]
        wavelength = d["wavelength"]
        template = d["template"]
        depths = d["depths"]
        avg_avg_depth = d["avg_avg_depth"]
        std_avg_depth = d["std_avg_depth"]
        avg_rms_depth = d["avg_rms_depth"]
        std_rms_depth = d["std_rms_depth"]

        # get the errors
        err_avg_depth = std_avg_depth ./ sqrt(Nloop)
        err_rms_depth = std_rms_depth ./ sqrt(Nloop)

        # figure out ylims
        ylims = (round(minimum(avg_rms_depth .- std_rms_depth) - 0.1, sigdigits=2),
                 round(maximum(avg_rms_depth .+ std_rms_depth) + 0.1, sigdigits=2))

        # create figure objects
        fig, ax1 = plt.subplots()

        # plot the results
        fig, ax1 = plt.subplots()
        ax1.errorbar(depths, avg_rms_depth, yerr=err_rms_depth, capsize=3.0, color="black", fmt=".")
        ax1.fill_between(depths, avg_rms_depth .- std_rms_depth, avg_rms_depth .+ std_rms_depth, color="tab:blue", alpha=0.3)
        ax1.fill_betweenx(range(0.0, ylims[2], length=5), zeros(5), zeros(5) .+ 0.25, hatch="/", fc="black", ec="white", alpha=0.15, zorder=0)

        # set labels, etc.
        txt = ("\${\\rm " * replace(template, "_" => "\\ ") * "}\$")
        ax1.set_title(txt)
        ax1.set_xlabel(L"{\rm Line\ Depth}")
        ax1.set_ylabel(L"{\rm RMS}_{\rm RV}\ {\rm (m s}^{-1})")
        ax1.set_xlim(0.0, 1.0)
        ax1.set_ylim(ylims...)


        # save and close the fig
        fig.savefig(plotdir * template * "_depth.pdf")
        plt.clf(); plt.close()
        println(">>> Figure written to: " * plotdir * template * "_depth.pdf")
    end
end


