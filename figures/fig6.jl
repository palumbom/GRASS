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
    # set up stuff for lines
    lines = [5434.5]
    depths = [0.8]
    resolution = 700000.0
    top = NaN
    contiguous_only=false

    # get angles
    n_inc = 20
    incls = (π/180.0) .* range(0.0, 90.0, length=n_inc)
    ycomp = sin.(incls)
    zcomp = cos.(incls)

    # make pole vectors
    poles = Array{Tuple{Float64, Float64, Float64}, 1}(undef, n_inc)
    for i in eachindex(incls)
        pole = (0.0, ycomp[i], zcomp[i])
        poles[i] = pole
    end

    # allocate shared sharred arrays
    avg_avg_inc = SharedArray{Float64}(n_inc)
    std_avg_inc = SharedArray{Float64}(n_inc)
    avg_rms_inc = SharedArray{Float64}(n_inc)
    std_rms_inc = SharedArray{Float64}(n_inc)

    # loop over inclinations
    @sync @distributed for i in eachindex(poles)
        println("running pole " * string(i) * " of " * string(n_inc))

         # create spec and disk params instances
        spec = SpecParams(lines=lines, depths=depths, resolution=resolution,
                          extrapolate=true, contiguous_only=contiguous_only)
        disk = DiskParams(N=N, Nt=Nt, pole=poles[i])

        # synthesize spectra, get velocities and stats
        avg_avg1, std_avg1, avg_rms1, std_rms1 = spec_loop(spec, disk, Nloop, top=top)
        avg_avg_inc[i] = avg_avg1
        std_avg_inc[i] = std_avg1
        avg_rms_inc[i] = avg_rms1
        std_rms_inc[i] = std_rms1
    end

    # make data frame
    df = DataFrame()
    df[!,:inc] = incls
    df[!,:avg_avg_inc] = avg_avg_inc
    df[!,:std_avg_inc] = std_avg_inc
    df[!,:avg_rms_inc] = avg_rms_inc
    df[!,:std_rms_inc] = std_rms_inc

    # write results to CSV
    fname = datadir * "inclination_" *  string(N) * ".csv"
    CSV.write(fname, df)
    return nothing
end

# run the simulation
if run
    main()
end

# plotting code block
if plot
    # plotting imports
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

    # read in the data
    fname = datadir * "inclination_" *  string(N) * ".csv"
    df = CSV.read(fname, DataFrame)

    # assign to variable names
    ang = df.inc .* (180.0/π)
    avg_avg_inc = df.avg_avg_inc
    std_avg_inc = df.std_avg_inc
    avg_rms_inc = df.avg_rms_inc
    std_rms_inc = df.std_rms_inc

    # get the errors
    err_avg_inc = std_avg_inc ./ sqrt(Nloop)
    err_rms_inc = std_rms_inc ./ sqrt(Nloop)

    # plot the results
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.errorbar(ang, avg_rms_inc, yerr=err_rms_inc, capsize=3.0, color="black", fmt=".")
    ax1.fill_between(ang, avg_rms_inc .- std_rms_inc, avg_rms_inc .+ std_rms_inc, color="tab:blue", alpha=0.3)

    # annotate the axes
    arrowprops = Dict("facecolor"=>"black", "shrink"=>0.05, "width"=>2.0,"headwidth"=>8.0)
    ax1.annotate(L"\textnormal{Pole-on}", xy=(69.8, 0.39), xytext=(0.0,0.388), arrowprops=arrowprops)
    ax1.annotate(L"\textnormal{Equator-on}", xy=(70.0, 0.388))

    # set labels, etc.
    ax1.set_xlabel(L"{\rm Inclination\ (deg)}")
    ax1.set_ylabel(L"{\rm RMS}_{\rm RV}\ {\rm (m s}^{-1})")
    ax1.set_xticks(range(0, 90, length=10))
    ax1.set_ylim(0.375, 0.675)

    # save the figure
    fig.savefig(plotdir * "fig6.pdf")
    plt.clf(); plt.close()
    println(">>> Figure written to: " * plotdir * "fig6.pdf")
end

