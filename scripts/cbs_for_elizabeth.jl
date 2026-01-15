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
lines_to_do = ["FeI_5434"]

# marker/colors for directions
cs = ["k", "tab:blue", "tab:orange", "tab:green", "tab:pink"]
ms = [".", "^", "v", "<", ">"]
lbls = ["Center", "North", "South", "East", "West"]

# Define the function
function model(μ, CB1, CB2, CB3)
    CB1 * exp(μ)^2.0 + CB2 * exp(μ) + CB3
end

# Define all coefficient sets with labels
datasets = [
    ("Fe I 5434 Å", [0.48680233581206817,-2.360771253495881,2.168676429569747]) #(CB + MF)
]

μ_vals = range(0, 1, length=1000)

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
            ax1.scatter([mu_val], [cbs], color=cs[ax_idx + 1], marker=ms[ax_idx+1], label=label) 
            seen[ax_idx+1] = true
        else
            ax1.scatter([mu_val], [cbs], color=cs[ax_idx + 1], marker=ms[ax_idx+1])
        end
    end

    # take the averages
    avgs ./= cnts
    ax1.scatter(mus_plot, avgs, c="k", marker="x", s=50, label="μ average")

    ax1.set_xlabel("μ", fontsize=12)
    ax1.set_ylabel("Convective Blueshift (m/s)", fontsize=12)
    # ax1.set_xlim(reverse(ax1.get_xlim())...)
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

    label, (CB1, CB2, CB3) = datasets[1]
    y_solid = [model(μ, CB1, CB2, CB3) for μ in μ_vals]
    ax1.plot(μ_vals, y_solid .* 1000 .+ 500, color = "k")
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(line_names[i] * "_cbs.pdf")
    # plt.clf(); plt.close()
end

