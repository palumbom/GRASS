# using Pkg; Pkg.activate(".")
using GRASS
using Printf
using Revise
using CSV, HDF5
using DataFrames
using Statistics
using EchelleCCFs
using ProgressMeter
using BenchmarkTools

# # plotting
# using LaTeXStrings
# import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
# mpl.style.use(GRASS.moddir * "fig.mplstyle")
# mpl.use("Qt5Agg")

# output directory
outdir = abspath("/mnt/ceph/users/mpalumbo/data_for_eduardo/")
!isdir(outdir) && mkdir(outdir)

# get the idx to run
template_idx = tryparse(Int, ARGS[1])

# get data for the template lines in GRASS library
lp = GRASS.LineProperties()
wavelength = lp.λrest
depth = lp.depth
dfile = lp.file
lname = GRASS.get_name(lp)

# set up parameters for synthetic spectrum
Nt = 1920 # number of 15-second time steps in simulation
disk = DiskParams(Nt=Nt)

let i = template_idx
    # parameters for synthetic spectrum
    lines = [wavelength[i]]
    depths = [depth[i]]
    templates = [dfile[i]]
    variability = trues(length(lines)) # control whether lines dance
    resolution = 7e5 # don't set resolution here, rather convolve it down later

    # spec object
    spec = SpecParams(lines=lines, depths=depths, variability=variability,
                      templates=templates, resolution=resolution, 
                      oversampling=4.0)

    # set mu bins
    μ_bins = range(0.15, 0.95, step=0.1)

    # synthesize
    wavs, flux = GRASS.synthesize_spectra_resolved(μ_bins, spec, disk, verbose=true, use_gpu=true)

    # write noiseless, full res spectra to disk
    fname = joinpath(outdir, lname[i] * "_noiseless.h5")
    h5open(fname, "w") do file
        write(file, "wavs", wavs, "flux", flux) 
    end

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

nbins = range(3, 10, step=1)
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
plt.ylabel(L"{\rm RMS\ RV\ of\ disk\ center\ bin\ (m/s)}")
plt.show()