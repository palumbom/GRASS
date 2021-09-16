# import stuff
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Statistics
@everywhere using GRASS
@everywhere using SharedArrays
@everywhere using EchelleCCFs
using CSV
using JLD2
using LsqFit
using FileIO
using DataFrames
using LaTeXStrings
using HypothesisTests

# define rms loop function
include(GRASS.moddir * "figures/fig_functions.jl")

# some global stuff
const N = 132
const Nloop = 1200

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

function main()
    # set observing parameters
    N_obs = range(2, 4, step=1)  # number of observations
    exp_time = 600.0             # ~2 p-mode periods
    snr = 5000.0                 # SNR per res element to get 50 cm/s precision
    new_res = 1.17e5

    # dead times
    short_dead = 30.0            # representative for NEID read-out time
    large_dead = 3600.0          # ~60 mins

    # set up parameters for spectrum
    lines = [5434.5]
    depths = [0.8]
    resolution = 7e5

    # allocate shared arrays for short off-target time
    avg_shortdead = SharedArray{Float64}(length(N_obs), Nloop)
    rms_shortdead = SharedArray{Float64}(length(N_obs), Nloop)

    # allocate shared arrays for large off-target time
    avg_largedead = SharedArray{Float64}(length(N_obs), Nloop)
    rms_largedead = SharedArray{Float64}(length(N_obs), Nloop)

    # create spec params instance
    spec1 = SpecParams(lines=lines, depths=depths, resolution=resolution)

    # loop over number of observations
    for i in eachindex(N_obs)
        # create observation plan instances
        obs_shortdead = GRASS.ObservationPlan(N_obs=N_obs[i],
                                              obs_per_night=N_obs[i],
                                              time_per_obs=exp_time,
                                              dead_time=short_dead)
        obs_largedead = GRASS.ObservationPlan(N_obs=N_obs[i],
                                              obs_per_night=N_obs[i],
                                              time_per_obs=exp_time,
                                              dead_time=large_dead)

        # compute vels repeatedly to get good stats
        @sync @distributed for j in 1:Nloop
            # first simulate short off-target time
            rvs1, sigs1 = GRASS.simulate_observations(obs_shortdead, spec1, N=N, snr=snr, new_res=new_res)
            avg_shortdead[i, j] = mean(rvs1)
            rms_shortdead[i, j] = calc_rms(rvs1)


            # then simulate large off-target time
            rvs2, sigs2 = GRASS.simulate_observations(obs_largedead, spec1, N=N, snr=snr, new_res=new_res)
            avg_largedead[i, j] = mean(rvs2)
            rms_largedead[i, j] = calc_rms(rvs2)
        end
    end

    # convert to plain arrays
    avg_shortdead = convert(Array{Float64}, avg_shortdead)
    rms_shortdead = convert(Array{Float64}, rms_shortdead)
    avg_largedead = convert(Array{Float64}, avg_largedead)
    rms_largedead = convert(Array{Float64}, rms_largedead)

    # save the output
    outfile = datadir * "observe_" * string(N) * "_loop_" * string(Nloop) * ".jld2"
    save(outfile,
         "avg_shortdead", avg_shortdead,
         "rms_shortdead", rms_shortdead,
         "avg_largedead", avg_largedead,
         "rms_largedead", rms_largedead,
         "exp_time", exp_time,
         "depths", depths,
         "N_obs", N_obs,
         "snr", snr)
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
    file = datadir * "observe_" * string(N) * "_loop_" * string(Nloop) * ".jld2"
    d = load(file)
    avg_shortdead = d["avg_shortdead"]
    rms_shortdead = d["rms_shortdead"]
    avg_largedead = d["avg_largedead"]
    rms_largedead = d["rms_largedead"]
    exp_time = d["exp_time"]
    depths = d["depths"]
    N_obs = d["N_obs"]
    snr = d["snr"]

    # initialize plotting objects
    fig = plt.figure(figsize=(10.5, 6.75))
    ax1 = plt.subplot2grid((2,4),(0,0), colspan=2)
    ax2 = plt.subplot2grid((2,4),(0,2), colspan=2)
    ax3 = plt.subplot2grid((2,4),(1,1), colspan=2)
    axs = [ax1, ax2, ax3]

    # make arrays for limits and loop over plotting windows
    xlims = []
    ylims = []
    for i in eachindex(axs)
        # histogram for short dead time
        axs[i].hist(avg_shortdead[i,:], density=true, histtype="stepfilled", fc="k", ec="k", lw=2, alpha=0.5, label="Short wait")
        axs[i].hist(avg_shortdead[i,:], density=true, histtype="step", ec="k", lw=1.0, alpha=0.8)

        # histogram for long wait time
        axs[i].hist(avg_largedead[i,:], density=true, histtype="stepfilled", fc="tab:blue", ec="tab:blue", lw=1.5, alpha=0.5, label="Long wait")
        axs[i].hist(avg_largedead[i,:], density=true, histtype="step", ec="tab:blue", lw=1.5, alpha=0.8)

        # get axis limits
        append!(xlims, [axs[i].get_xlim()])
        append!(ylims, [axs[i].get_ylim()])

        # report two-sample KS test
        @show(ApproximateTwoSampleKSTest(avg_shortdead[i,:], avg_largedead[i,:]))
        println()
    end

    # set the axes limits
    bbox = Dict("ec"=>"black", "fc"=>"white")
    xy = (maximum(last.(xlims)) - 0.25, maximum(last.(ylims)) - 0.4)
    for i in eachindex(axs)
        # annotate axes and set labels, etc.
        axs[i].annotate(L"N_{\rm obs} =\ " * latexstring(N_obs[i]), xy=xy, bbox=bbox)
        axs[i].set_xlabel(L"{\rm Mean\ velocity\ (m s}^{-1})")
        axs[i].set_ylabel(L"{\rm Probability\ density}")
        axs[i].set_xlim(round(minimum(first.(xlims)), digits=1),
                        round(maximum(last.(xlims)), digits=1))
        axs[i].set_ylim(minimum(first.(ylims)), maximum(last.(ylims)))
        axs[i].xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    end

    # set legend in one panel and write figure
    axs[1].legend(loc="upper left")
    fig.savefig(plotdir * "fig7.pdf")
    plt.clf(); plt.close()
    println(">>> Figure written to: " * plotdir * "fig7.pdf")
end
