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
include(GRASS.moddir * "figures/rms_loop.jl")

# some global stuff
const N = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]
const Nt = 100
const Nloop = 12

function resolution()
    # set up parameters for lines
    lines = [5434.5]
    depths = [0.8]
    res = 700000.0
    top = NaN
    contiguous_only=true
    spec = SpecParams(lines=lines, depths=depths, resolution=res, extrapolate=true,
                      contiguous_only=contiguous_only)

    # allocate shared arrays
    rms_res = SharedArray{Float64}(length(N))
    std_res = SharedArray{Float64}(length(N))

    # calculate
    @sync @distributed for i in 1:length(N)
    	println("running resolution N = " * string(N[i]))
        disk = DiskParams(N=N[i], Nt=Nt)
    	rms1, std1 = rms_loop(spec, disk, Nloop, top=top)
    	rms_res[i] = rms1
    	std_res[i] = std1
    end

    # make data frame
    df = DataFrame()
    df[!,:res] = N
    df[!,:rms_res] = rms_res
    df[!,:std_res] = std_res

    # write to CSV
    fname = SS.moddir * "scripts/out/rms_vs_res.csv"
    CSV.write(fname, df)
    return nothing
end

# run only if on ACI, else just plot the saved results
if SS.moddir == "/storage/work/m/mlp95/SynthSpectra/"
    resolution()
    plot_resolution = false
else
    plot_resolution = true
end

# plotting code block
if plot_resolution
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use("my.mplstyle")

    # set Desktop directory
    outdir = "/Users/michael/Desktop/"
    if !isdir(outdir)
        outdir = "/Users/mlp95/Desktop/"
    end

    # read in the data
    fname = SS.moddir * "scripts/out/rms_vs_res.csv"
    df = CSV.read(fname, DataFrame)
    res = df.res
    rms_res = df.rms_res
    std_res = df.std_res

    # fit the data
    @. power_law(x, p) = p[1] * x^(-p[2])
    fit = curve_fit(power_law, res, rms_res, [1.0, 1.0])
    res_fit = range(6, 2400, length=1000)
    println("Best fit power law index = " * string(fit.param[2]))

    x = string(round(fit.param[2], sigdigits=3))

    # plot it
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log", base=10)
    ax1.errorbar(res, rms_res, yerr=std_res, capsize=3.0, color="black", fmt=".")
    ax1.plot(res_fit, power_law(res_fit, fit.param), "k--", alpha=0.4,
             label = L"{\rm Power\ law\ index\ } \approx\ %$x ")

    xrng = ax1.get_xlim()

    ax1.fill_between(range(xrng[1], xrng[2], length=2), repeat([0.319-0.09], 2), repeat([0.319+0.09], 2),
                     alpha=0.5, color="tab:blue", label=L"{\rm Elsworth\ et\ al.\ (1994)}")
    ax1.fill_between(range(xrng[1], xrng[2], length=2), repeat([0.461-0.1], 2), repeat([0.461+0.1], 2),
                     alpha=0.5, color="tab:green", label=L"{\rm Palle\ et\ al.\ (1999)}")
    ax1.set_xlim(xrng...)
    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.set_xlabel(L"N")
    ax1.set_ylabel(L"{\rm RMS\ RV\ (m s}^{-1})")
    ax1.legend()
    fig.savefig(outdir * "res.pdf")
    plt.clf(); plt.close()
end
