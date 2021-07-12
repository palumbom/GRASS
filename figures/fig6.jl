# import stuff
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Statistics
@everywhere using SynthSpectra
@everywhere SS=SynthSpectra
@everywhere using SharedArrays
@everywhere using EchelleCCFs
using CSV
using DataFrames
using LaTeXStrings

# define rms loop
include(string(@__DIR__) * "/rms_loop.jl")

# some global stuff
const N = 256
const Nt = 200
const Nloop = 100

function inclination()
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
    prms = SharedArray{Float64}(n_inc)
    pstd = SharedArray{Float64}(n_inc)

    # calculate
    @sync @distributed for i in eachindex(poles)
        println("running pole " * string(i) * " of " * string(n_inc))
        spec = SpecParams(lines=lines, depths=depths, resolution=resolution,
                          extrapolate=true, contiguous_only=contiguous_only)
        disk = DiskParams(N=N, Nt=Nt, pole=poles[i])
        rms1, std1 = rms_loop(spec, disk, Nloop, top=top)
        prms[i] = rms1
        pstd[i] = std1
    end

    # make data frame
    df = DataFrame()
    df[!,:inc] = incls
    df[!,:rms] = prms
    df[!,:std] = pstd

    # write results to CSV
    fname = SS.moddir * "scripts/out/inclination_" *  string(N) * ".csv"
    CSV.write(fname, df)
    return nothing
end

# run only if on ACI, else just plot the saved results
if SS.moddir == "/storage/work/m/mlp95/SynthSpectra/"
    inclination()
    plot_inclination = false
else
    plot_inclination = true
end

# plotting code block
if plot_inclination
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use("my.mplstyle")

    # set Desktop directory
    outdir = "/Users/michael/Desktop/"
    if !isdir(outdir)
        outdir = "/Users/mlp95/Desktop/"
    end

    # read in the data
    fname = SS.moddir * "scripts/out/inclination_" *  string(N) * ".csv"
    dat = CSV.read(fname, DataFrame)
    ang = dat.inc .* (180.0/π)
    rms = dat.rms
    tstd = dat.std

    arrowprops=Dict("facecolor"=>"black", "shrink"=>0.05, "width"=>2.0,"headwidth"=>8.0)

    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.errorbar(ang, rms, yerr=tstd, capsize=3.0, color="black", fmt=".")
    ax1.set_xlabel(L"{\rm Inclination\ (deg)}")
    ax1.set_ylabel(L"{\rm RMS\ RV\ (m s}^{-1})")
    ax1.set_xticks(range(0, 90, length=10))
    ax1.set_ylim(0.2,0.345)
    ax1.annotate("Pole-on", xy=(77.5, 0.207), xytext=(0.5,0.205), arrowprops=arrowprops)
    ax1.annotate("Equator-on", xy=(78.0, 0.205))
    fig.savefig(outdir * "inc_rms.pdf")
    plt.clf(); plt.close()
end

