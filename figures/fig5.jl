# import stuff
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Statistics
@everywhere using GRASS
@everywhere using SharedArrays
@everywhere using EchelleCCFs
using CSV
using DataFrames
using LaTeXStrings
using LsqFit

# define rms loop function
include(GRASS.moddir * "figures/fig_functions.jl")

# some global stuff
const N = 132
const Nt = 200
const Nloop = 200

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

function main()
    # set up parameters for lines
    lines = [5434.5]
    depths = range(0.05, stop=0.95, step=0.05)
    resolution=700000.0
    top = NaN
    contiguous_only=false

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
        spec = SpecParams(lines=lines, depths=[depths[i]], resolution=resolution,
                          extrapolate=true, contiguous_only=contiguous_only)

        # synthesize spectra, get velocities and stats
        avg_avg1, std_avg1, avg_rms1, std_rms1 = spec_loop(spec, disk, Nloop, top=top)
        avg_avg_depth[i] = avg_avg1
        std_avg_depth[i] = std_avg1
        avg_rms_depth[i] = avg_rms1
        std_rms_depth[i] = std_rms1
    end

    # make data frame
    df = DataFrame()
    df[!,:depths] = depths
    df[!,:avg_avg_depth] = avg_avg_depth
    df[!,:std_avg_depth] = std_avg_depth
    df[!,:avg_rms_depth] = avg_rms_depth
    df[!,:std_rms_depth] = std_rms_depth

    # write to CSV
    fname = datadir * "rms_vs_depth_" * string(N) * ".csv"
    CSV.write(fname, df)
    return nothing
end

# run the simulation
if run
    main()
end

if plot
    # plotting imports
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

    # read in the data
    fname = datadir * "rms_vs_depth_" * string(N) * ".csv"
    df = CSV.read(fname, DataFrame)

    # assign to variable names
    depths = df.depths
    avg_avg_depth = df.avg_avg_depth
    std_avg_depth = df.std_avg_depth
    avg_rms_depth = df.avg_rms_depth
    std_rms_depth = df.std_rms_depth

    # get the errors
    err_avg_depth = std_avg_depth ./ sqrt(Nloop)
    err_rms_depth = std_rms_depth ./ sqrt(Nloop)

    # plot the results
    fig, ax1 = plt.subplots()
    ax1.errorbar(depths, avg_rms_depth, yerr=err_rms_depth, capsize=3.0, color="black", fmt=".")
    ax1.fill_between(depths, avg_rms_depth .- std_rms_depth, avg_rms_depth .+ std_rms_depth, color="tab:blue", alpha=0.3)
    ax1.fill_betweenx(range(0.0, 1.0, length=5), zeros(5), zeros(5) .+ 0.2, hatch="/", fc="black", ec="white", alpha=0.15, zorder=0)

    # set labels, etc.
    ax1.set_xlabel(L"{\rm Line\ Depth}")
    ax1.set_ylabel(L"{\rm RMS}_{\rm RV}\ {\rm (m s}^{-1})")
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(0.375, 0.675)

    # annotate the axes and save the figure
    arrowprops = Dict("facecolor"=>"black", "shrink"=>0.05, "width"=>2.0,"headwidth"=>8.0)
    ax1.annotate(L"{\rm Shallow}", xy=(0.85,0.39), xytext=(0.05,0.388), arrowprops=arrowprops)
    ax1.annotate(L"{\rm Deep}", xy=(0.86, 0.388))
    fig.savefig(plotdir * "fig5.pdf")
    plt.clf(); plt.close()
    println(">>> Figure written to: " * plotdir * "fig5.pdf")
end


