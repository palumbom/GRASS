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
using JLD2
using FileIO
using LaTeXStrings

# define rms loop function
include(GRASS.moddir * "figures/fig_functions.jl")

# some global stuff
const N = 256
const Nloop = 800

# set plotting boolean
plot = true

# check directories
grassdir, plotdir, datadir = check_plot_dirs()

function observe()
    # set observing parameters
    N_obs = range(2, 4, step=1)  # ~20 minutes continuous at most
    exp_time = 600.0             # ~2 p-mode periods
    snr = 3500.0                 # SNR per res element to get 50 cm/s precision
    new_res = 1.17e5

    # dead times
    short_dead = 30.0            # representative for NEID read-out time
    large_dead = 3600.0          # ~60 mins

    # set up parameters for spectrum
    lines = [5434.5]
    depths = [0.8]
    resolution = 7e5

    # allocate shared arrays
    avg_shortdead = SharedArray{Float64}(length(N_obs))
    avg_largedead = SharedArray{Float64}(length(N_obs))
    rms_shortdead = SharedArray{Float64}(length(N_obs))
    rms_largedead = SharedArray{Float64}(length(N_obs))

    # create spec params instance
    spec1 = SpecParams(lines=lines, depths=depths, resolution=resolution)

    # now loop over and create observation plans
    @sync @distributed for i in eachindex(N_obs)
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
            vels1 = GRASS.simulate_observations(obs_shortdead, spec1, N=N, snr=snr, new_res=new_res)
            vels2 = GRASS.simulate_observations(obs_largedead, spec1, N=N, snr=snr, new_res=new_res)
            avg_shortdead[i] += mean(vels1)
            avg_largedead[i] += mean(vels2)
            rms_shortdead[i] += calc_rms(vels1)
            rms_largedead[i] += calc_rms(vels2)
        end
    end

    # compute the mean values
    avg_shortdead[i] ./= Nloop
    avg_largedead[i] ./= Nloop
    rms_shortdead[i] ./= Nloop
    rms_largedead[i] ./= Nloop

    # make into regular arrays
    avg_shortdead = convert(Array{Float64}, avg_shortdead)
    avg_largedead = convert(Array{Float64}, avg_largedead)
    rms_shortdead = convert(Array{Float64}, rms_shortdead)
    rms_largedead = convert(Array{Float64}, rms_largedead)

    # save the output
    outfile = datadir * "observe_" * string(N) * "_loop_" * string(Nloop) * ".jld2"
    save(outfile,
         "avg_shortdead", avg_shortdead,
         "avg_largedead", avg_largedead,
         "rms_shortdead", rms_shortdead,
         "rms_largedead", rms_largedead,
         "exp_time", exp_time,
         "depths", depths,
         "N_obs", N_obs,
         "snr", snr)
    return nothing
end

# run the simulation
observe()

# plotting code block
if plot
    # plotting imports
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

    # other imports for plotting
    using LsqFit
    using HypothesisTests

    # read in the data
    outfile = datadir * "observe_" * string(N) * "_loop_" * string(Nloop) * ".jld2"
    d = load(datadir)
    avg_shortdead = d["avg_shortdead"]
    avg_largedead = d["avg_largedead"]
    # rms_shortdead = d["rms_shortdead"]
    # rms_largedead = d["rms_largedead"]
    exp_time = d["exp_time"]
    depths = d["depths"]
    N_obs = d["N_obs"]
    snr = d["snr"]

    # for i in eachindex(depths)
    let i = 4
        fig = plt.figure(figsize=(9.75,6.75))
        ax1 = plt.subplot2grid((2,4),(0,0), colspan=2)
        ax2 = plt.subplot2grid((2,4),(0,2), colspan=2)
        ax3 = plt.subplot2grid((2,4),(1,1), colspan=2)
        axs = [ax1, ax2, ax3]

        xlims = []
        for j in eachindex(axs)
            axs[j].hist(avg_shortdead[i,j,:], density=true, histtype="stepfilled", fc="k", ec="k", lw=2, alpha=0.5, label="Short wait")
            axs[j].hist(avg_shortdead[i,j,:], density=true, histtype="step", ec="k", lw=1.0, alpha=0.8)

            axs[j].hist(avg_largedead[i,j,:], density=true, histtype="stepfilled", fc="tab:blue", ec="tab:blue", lw=1.5, alpha=0.5, label="large wait")
            axs[j].hist(avg_largedead[i,j,:], density=true, histtype="step", ec="tab:blue", lw=1.5, alpha=0.8)

            axs[j].annotate(L"N_{\rm obs} =\ " * latexstring(N_obs[j]), xy=(-204.15, 2.5), backgroundcolor="white")
            axs[j].set_xlabel(L"{\rm Measured\ velocity\ (m s}^{-1})")
            axs[j].set_ylabel(L"{\rm Probability\ density}")
            axs[j].set_ylim(0,2.75)
            append!(xlims, [axs[j].get_xlim()])

            # report two-sample KS test
            @show(ApproximateTwoSampleKSTest(avg_shortdead[i,j,:], avg_largedead[i,j,:]))
        end

        # set the xlimits and other nice things
        for j in eachindex(axs)
            axs[j].set_xlim(minimum(first.(xlims)), maximum(last.(xlims)))
        end
        axs[1].legend(loc="upper left")
        # axs[2].yaxis.set_label_position("right")
        # axs[2].yaxis.tick_right()
        fig.savefig(plotdir * "fig7.pdf")
        plt.clf(); plt.close()
        println(">>> Figure written to: " * plotdir * "fig7.pdf")
    end
end
