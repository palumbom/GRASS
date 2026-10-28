using CSV
using CUDA
using GRASS
using Printf
using Revise
using PyCall
using DataFrames
using Statistics
using EchelleCCFs
using Polynomials
using BenchmarkTools
using HypothesisTests

# plotting
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
pe = pyimport("matplotlib.patheffects")

# optional: define colors manually if needed (since mplstyle is gone)
color = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

# get data
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
files = GRASS.get_file(lp)
line_names = GRASS.get_name(lp)

# get lines to do 
lines_to_do = ["FeI_6173"]
a = 1.4097714740057676
b = -2.0650721844470294

df = CSV.read("data/disco_Fe6173_cb.csv", DataFrame)
DISCOCB = df.DISCOCB
PolynomialCB = df.PolynomialCB
MeasuredCB = df.MeasuredCB
mus_disco = df.Mu
# marker/colors for directions
cs = ["k", "tab:blue", "tab:orange", "tab:green", "tab:pink"]
ms = [".", ".", ".", ".", "."]
lbls = ["Center", "North", "South", "East", "West"]

# loop over files
for i in eachindex(files)
    !(line_names[i] in lines_to_do) && continue
    @show i

    # read in the data 
    soldata = GRASS.SolarData(fname=files[i])

    # make arrays 
    mus = unique(soldata.mu)
    mus_plot = zeros(length(mus))
    avgs = zeros(length(mus))
    cnts = zeros(length(mus))
    seen = falses(length(ms))

    # make a figure
    fig, ax1 = plt.subplots()

    # loop over keys 
    for (j,k) in enumerate(keys(soldata.cbs))
        # parse out axis and mu
        ax_idx = GRASS.parse_ax_string(string(k[1]))
        mu_val = GRASS.parse_mu_string(string(k[2]))
        
        cbs = soldata.cbs[k] * GRASS.c_ms

        # get index for the mu 
        mu_idx = findfirst(mu_val .== mus)
        avgs[mu_idx] += cbs
        cnts[mu_idx] += 1
        mus_plot[mu_idx] = mu_val

        if !seen[ax_idx+1]
            label = lbls[ax_idx + 1]
            ax1.errorbar(
                                [mu_val], [cbs];
                                yerr=[5.0],
                                fmt=ms[ax_idx+1],
                                color=cs[ax_idx+1],
                                capsize=3,
                                label=label
                            ) 
            seen[ax_idx+1] = true
        else
            ax1.errorbar(
                            [mu_val], [cbs];
                            yerr=[5.0],
                            fmt=ms[ax_idx+1],
                            color=cs[ax_idx+1],
                            capsize=3
                        )
        end
    end

    # generate smooth curve
    xfit = range(0, 1, length=500)

    yfit = a .* xfit.^2 .+ b .* xfit

    ax1.plot(xfit, (yfit * 1000) .- mean(yfit * 1000), color = "k", lw=2)

    ax1.plot(mus_disco, PolynomialCB, color = "g", lw=2, label = "Polynomial CB")
    ax1.plot(mus_disco, PolynomialCB .- mean(PolynomialCB), color = "g", ls = "--", lw=2)
    ax1.plot(mus_disco, MeasuredCB, color = "b", lw=2, label = "Measured CB")
    ax1.plot(mus_disco, DISCOCB, color = "pink", lw=2, label = "DISCO CB")

    ax1.set_xlabel("μ", fontsize=12)
    ax1.set_ylabel("Convective Blueshift (m/s)", fontsize=12)
    ax1.invert_xaxis()
    ax1.legend(ncols=2, fontsize=12)

    # generate title
    title = replace(line_names[i], "_" => " ")
    idx = findfirst('I', title)
    if isnothing(idx)
        line_title = title
    else
        line_title = title[1:idx-1] * " " * title[idx:end] * " Å"
    end
    ax1.set_title(line_title, fontsize=12)
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax1.set_xlim(1.02, 0)
    plt.savefig(line_names[i] * "_cbs.png", dpi = 600, bbox_inches="tight")
end

