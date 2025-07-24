# using Pkg; Pkg.activate(".")
using GRASS
using Printf
using Revise
using Statistics
using EchelleCCFs
using ProgressMeter
using BenchmarkTools

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")
mpl.use("Qt5Agg")

# set up paramaters for spectrum
N = 197
Nt = 500

lines = [5432.0]#, 5432.9]
depths = [0.6]#, 0.75]
templates = ["FeI_5432"]#, "FeI_6173"]
variability = trues(length(lines))
blueshifts = zeros(length(lines))
resolution = 7e5
seed_rng = false

disk = DiskParams(N=N, Nt=Nt, inclination=90.0, Nsubgrid=40)#, u1=0.3033, u2=0.3133, vsini=33950.0, A=0.7, B=0.0, C=0.0)
spec = SpecParams(lines=lines, depths=depths, variability=variability,
                   blueshifts=blueshifts, templates=templates,
                   resolution=resolution, buffer=0.5)

# set mu bins
μ_bins = range(0.15, 0.95, step=0.1)

# do the synthesis
println(">>> Synthesizing on GPU...")
tstart = time()
lambdas_gpu, outspec_gpu = GRASS.synthesize_spectra_resolved(μ_bins, spec, disk, seed_rng=seed_rng, verbose=true, use_gpu=true)
tstop = time()
@printf(">>> Synthesis time --> %.3f seconds \n", tstop - tstart)

# for μ_idx in eachindex(μ_bins)
#     # plt.plot(lambdas_gpu, outspec_gpu[:,μ_idx,1] ./ maximum(outspec_gpu[:,μ_idx,1]))
#     plt.plot(lambdas_gpu, outspec_gpu[:,μ_idx,1])
# end
# plt.plot(lambdas_gpu, sum(outspec_gpu[:,:,1], dims=2), c="k")
# plt.show()

# measure velocities
# rvs_out = zeros(Nt, length(μ_bins))
# for μ_idx in eachindex(μ_bins)
#     outspec_view = view(outspec_gpu, :, μ_idx, :)

#     v_grid_gpu, ccf_gpu = calc_ccf(lambdas_gpu, outspec_view, spec)
#     rvs_gpu, sigs_gpu = calc_rvs_from_ccf(v_grid_gpu, ccf_gpu)

#     rvs_out[:,μ_idx] .= rvs_gpu

#     plt.scatter(1:Nt, rvs_gpu)#.- mean(rvs_gpu))
# end
# plt.show()

nbins = range(3, 5, step=1)
ntrials = 10
@showprogress for nn in eachindex(nbins)
    rms_rv = zeros(ntrials)
    for jj in 1:ntrials
        # make bins
        μ_bins1 = range(0.1, 0.9, length=nbins[nn])

        # make spectra
        lambdas_gpu1, outspec_gpu1 = GRASS.synthesize_spectra_resolved(μ_bins1, spec, disk, seed_rng=seed_rng, verbose=false, use_gpu=true)

        μ_idx = lastindex(μ_bins1)
        outspec_view = view(outspec_gpu1, :, μ_idx, :)

        v_grid_gpu, ccf_gpu = calc_ccf(lambdas_gpu1, outspec_view, spec)
        rvs_gpu, sigs_gpu = calc_rvs_from_ccf(v_grid_gpu, ccf_gpu)

        rms_rv[jj] = calc_rms(rvs_gpu)
    end

    plt.errorbar([nbins[nn]], [mean(rms_rv)], yerr=[std(rms_rv)], c="k")
end
plt.xlabel(L"{\rm Number\ of\ \mu\ bins}")
plt.ylabel(L"{\rm RMS\ RV\ of\ disk-center\ bin\ (m/s)}")
plt.show()