# import stuff
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Statistics
@everywhere using CUDA
@everywhere using GRASS
@everywhere using SharedArrays
@everywhere using EchelleCCFs
using CSV
using Glob
using DataFrames
using LaTeXStrings
using LsqFit

# define rms loop function
include(GRASS.moddir * "figures1/fig_functions.jl")

# some global stuff
const N = 132
const Nt = 200
const Nloop = 100

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# decide whether to use gpu
use_gpu = CUDA.functional()

# lines to run
# airwavs = [5250.2084, 5432.9470, 5434.5232, 5576.0881, 6173.3344]
airwavs = [5432.9470, 5434.5232, 6173.3344]
# indirs = ["FeI_5250.2/", "FeI_5432/", "FeI_5434/", "FeI_5576/", "FeI_6173/"]
indirs = ["FeI_5432/", "FeI_5434/", "FeI_6173/"]

function single_line_variability(airwav, indir)
    # set up parameters for lines
    lines = [airwav]
    indirs = [GRASS.soldir * indir]
    depths = range(0.05, stop=0.95, step=0.05)
    resolution=700000.0
    top = NaN

    # allocate shared arrays
    avg_avg_depth = SharedArray{Float64}(length(depths))
    std_avg_depth = SharedArray{Float64}(length(depths))
    avg_rms_depth = SharedArray{Float64}(length(depths))
    std_rms_depth = SharedArray{Float64}(length(depths))

    # loop over depths
    disk = DiskParams(N=N, Nt=Nt)
    @sync @distributed for i in eachindex(depths)
        println("running depth = " * string(depths[i]))

        # create spec instance
        spec = SpecParams(lines=lines, depths=[depths[i]], indirs=indirs, resolution=resolution)

        # synthesize spectra, get velocities and stats
        out = spec_loop(spec, disk, Nloop, use_gpu=use_gpu)
        avg_avg_depth[i] = out[1]
        std_avg_depth[i] = out[2]
        avg_rms_depth[i] = out[3]
        std_rms_depth[i] = out[4]
    end

    # make data frame
    df = DataFrame()
    df[!,:airwav] = repeat([airwav], length(depths))
    df[!,:depths] = depths
    df[!,:avg_avg_depth] = avg_avg_depth
    df[!,:std_avg_depth] = std_avg_depth
    df[!,:avg_rms_depth] = avg_rms_depth
    df[!,:std_rms_depth] = std_rms_depth

    # write to CSV
    fname = datadir * "rms_vs_depth_" * indir[1:end-1] * ".csv"
    CSV.write(fname, df)
    return nothing
end

# run the simulation
# if run
#     for l in eachindex(airwavs)
#         single_line_variability(airwavs[l], indirs[l])
#     end
# end

if plot
    # plotting imports
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

    # create data frame
    df = DataFrame(airwav=[], depths=[], avg_avg_depth=[], std_avg_depth=[], avg_rms_depth=[], std_rms_depth=[])

    # read in the data
    files = Glob.glob("rms_vs_depth_*.csv", datadir)
    for f in files
        if splitpath(f)[end] ==  "rms_vs_depth_FeI_5434.csv"
            continue
        elseif splitpath(f)[end] ==  "rms_vs_depth_132.csv"
            df_temp = CSV.read(f, DataFrame)
            df_temp[!, :airwav] = repeat([5434.5232], length(df_temp.depths))
            append!(df, df_temp)
            continue
        end
        append!(df, CSV.read(f, DataFrame))
    end

    # make sure its sorted on airwavs
    sort!(df, :airwav)

    # set color, label lists, etc.
    goodwavs = [5434.5232]#, 5432.9470, 6173.3344]
    geffs = [L"g_{\rm eff} = 0.00"]#, L"g_{\rm eff} = 0.50", L"g_{\rm eff} = 2.50"]
    colors = ["tab:blue"]#, "tab:green", "tab:orange"]

    # create figure objects
    fig, ax1 = plt.subplots()

    # loop over unique airwavs
    airwavs = unique(df.airwav)
    for i in eachindex(airwavs)
        if !(airwavs[i] in goodwavs)
            continue
        end

        # get index of goodwav
        j = findfirst(airwavs[i] .== goodwavs)

        # assign to variable names
        inds = df.airwav .== airwavs[i]
        depths = df.depths[inds]
        avg_avg_depth = df.avg_avg_depth[inds]
        std_avg_depth = df.std_avg_depth[inds]
        avg_rms_depth = df.avg_rms_depth[inds]
        std_rms_depth = df.std_rms_depth[inds]

        # get the errors
        err_avg_depth = std_avg_depth ./ sqrt(Nloop)
        err_rms_depth = std_rms_depth ./ sqrt(Nloop)

        # plot the results
        ax1.errorbar(depths, avg_rms_depth, yerr=err_rms_depth, capsize=3.0, color=colors[j], fmt=".", label=geffs[j])
        ax1.fill_between(depths, avg_rms_depth .- std_rms_depth, avg_rms_depth .+ std_rms_depth, color=colors[j], alpha=0.3)
    end

    # shade the low depth area
    ax1.fill_betweenx(range(0.0, 1.0, length=5), zeros(5), zeros(5) .+ 0.2, hatch="/", fc="black", ec="white", alpha=0.15, zorder=0)

    # set labels, etc.
    ax1.set_xlabel(L"{\rm Line\ Depth}")
    ax1.set_ylabel(L"{\rm RMS}_{\rm RV}\ {\rm (m s}^{-1})")
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(0.32, 0.9)

    # annotate the axes and save the figure
    arrowprops = Dict("facecolor"=>"black", "shrink"=>0.05, "width"=>2.0,"headwidth"=>8.0)
    ax1.annotate(" ", xy=(0.85,0.39), xytext=(0.2,0.39), arrowprops=arrowprops)
    ax1.annotate(L"{\rm Shallow}", xy=(0.05, 0.385))
    ax1.annotate(L"{\rm Deep}", xy=(0.86, 0.385))
    ax1.legend(loc="upper center", ncol=3, fontsize=12)
    fig.savefig(plotdir * "depths.pdf")
    plt.clf(); plt.close()
    println(">>> Figure written to: " * plotdir * "depths.pdf")
end


