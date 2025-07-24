# using Pkg; Pkg.activate(".")
using GRASS
using Printf
using Revise
using Statistics
using EchelleCCFs
using BenchmarkTools

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")
mpl.use("Qt5Agg")

# set up paramaters for spectrum
N = 197
Nt = 50

lines = [5432.0]#, 5432.9]
depths = [0.6]#, 0.75]
templates = ["FeI_5432"]#, "FeI_6173"]
variability = falses(length(lines))
blueshifts = zeros(length(lines))
resolution = 7e5
seed_rng = true

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

for μ_idx in eachindex(μ_bins)
    # plt.plot(lambdas_gpu, outspec_gpu[:,μ_idx,1] ./ maximum(outspec_gpu[:,μ_idx,1]))
    plt.plot(lambdas_gpu, outspec_gpu[:,μ_idx,1])
end
plt.plot(lambdas_gpu, sum(outspec_gpu[:,:,1], dims=2), c="k")
plt.show()

# measure velocities
rvs_out = zeros(Nt, length(μ_bins))
for μ_idx in eachindex(μ_bins)
    outspec_view = view(outspec_gpu, :, μ_idx, :)

    v_grid_gpu, ccf_gpu = calc_ccf(lambdas_gpu, outspec_view, spec)
    rvs_gpu, sigs_gpu = calc_rvs_from_ccf(v_grid_gpu, ccf_gpu)

    rvs_out[:,μ_idx] .= rvs_gpu
end