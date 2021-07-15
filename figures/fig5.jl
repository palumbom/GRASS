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

# plotting packages
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use("GRASS.moddir * "figures/fig.mplstyle"")

# define rms loop function
include(GRASS.moddir * "figures/rms_loop.jl")

# some global stuff
const N = 256
const Nt = 200
const Nloop = 100

function depth(;write=false)
    # set up parameters for lines
    lines = [5434.5]
    depths = range(0.05, stop=0.95, step=0.05)
    resolution=700000.0
    top = NaN
    contiguous_only=false

    # allocate shared arrays
    depth_rms_yp = SharedArray{Float64}(length(depths))
    depth_std_yp = SharedArray{Float64}(length(depths))

    # calculate
    disk = DiskParams(N=N, Nt=Nt)
    @sync @distributed for i in eachindex(depths)
        println("running depth = " * string(depths[i]))
        spec2 = SpecParams(lines=lines, depths=[depths[i]], resolution=resolution,
                           extrapolate=true, contiguous_only=contiguous_only)
        rms2, std2 = rms_loop(spec2, disk, Nloop, top=top)
        depth_rms_yp[i] = rms2
        depth_std_yp[i] = std2
    end

    # make data frame
    df = DataFrame()
    df[!,:depths] = depths
    df[!,:rms_yp] = depth_rms_yp
    df[!,:std_yp] = depth_std_yp

    # write to CSV
    if write
        fname = GRASS.moddir * "scripts/out/rms_vs_depth_" * string(N) * ".csv"
        CSV.write(fname, df)
    end
    return df
end

# run the simulation
df = depth(write=false)

# plot the results
arrowprops = Dict("facecolor"=>"black", "shrink"=>0.05, "width"=>2.0,"headwidth"=>8.0)
fig, ax1 = plt.subplots()
# ax1.errorbar(df.depths, df.rms_yp, yerr=df.std_yp, capsize=3.0, color="black", fmt=".")
ax1.plot(df.depths, df.rms_yp, color="black", fmt=".")
ax1.fill_betweenx(range(0.0, 1.0, length=5), zeros(5), zeros(5) .+ 0.2, color="k", alpha=0.25)
ax1.set_xlabel(L"{\rm Line\ Depth}")
ax1.set_ylabel(L"{\rm RMS\ RV\ (m s}^{-1})")
ax1.set_xlim(0.0, 1.0)
ax1.set_ylim(0.19, 0.36)
ax1.annotate("Shallow", xy=(0.85,0.202), xytext=(0.05,0.2), arrowprops=arrowprops)
ax1.annotate("Deep", xy=(0.86, 0.2))
fig.savefig("rms_vs_depth.pdf")
plt.clf(); plt.close()


