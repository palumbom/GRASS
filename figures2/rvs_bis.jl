using Pkg; Pkg.activate(".")
using CUDA
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

# set up paramaters for spectrum
N = 132
Nt = 1000
lines = [5434.5]#, 5436.6]
depths = [0.9]#, 0.65]
geffs = [0.0]#, 0.0]
templates = ["TiII_5381"]#, "FeI_5436.6"]
variability = repeat([true], length(lines))
resolution = 7e5
seed_rng = true

# FeI 5379
# Ti II 5381

disk = DiskParams(N=N, Nt=Nt)
spec1 = SpecParams(lines=lines, depths=depths, variability=variability,
                   geffs=geffs, templates=templates, resolution=resolution)
lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=false, verbose=true, use_gpu=true)

# loop over resolutions
oversampling = 3.0
resolutions = reverse([0.5e5, 1e5, 2e5, 7e5])
# resolutions = [2e5]
colors = ["k", "tab:blue", "tab:orange", "tab:green"]
# colors = ["tab:blue"]
for i in eachindex(resolutions)
    @show resolutions[i]
    if resolutions[i] == 7e5
        wavs_out = lambdas1
        flux_out = outspec1
    else
        # do an initial conv to get output size
        wavs_to_deg = view(lambdas1, :, 1)
        flux_to_deg = view(outspec1, :, 1)
        wavs_degd, flux_degd = GRASS.convolve_gauss(wavs_to_deg,
                                                    flux_to_deg,
                                                    new_res=resolutions[i],
                                                    oversampling=oversampling)

        # allocate memory
        wavs_out = zeros(size(wavs_degd, 1), size(outspec1, 2))
        flux_out = zeros(size(wavs_degd, 1), size(outspec1, 2))

        # loop over epochs and convolve
        for j in 1:size(outspec1,2)
            flux_to_deg = view(outspec1, :, j)
            wavs_degd, flux_degd = GRASS.convolve_gauss(wavs_to_deg,
                                                        flux_to_deg,
                                                        new_res=resolutions[i],
                                                        oversampling=oversampling)

            # copy to array
            wavs_out[:, j] .= wavs_degd
            flux_out[:, j] .= flux_degd
        end
    end

    v_grid, ccf1 = calc_ccf(wavs_out[:,1], flux_out, [5434.5], [0.9],
                            resolutions[i], normalize=true,
                            mask_type=EchelleCCFs.TopHatCCFMask)
    rvs1, sigs1 = calc_rvs_from_ccf(v_grid, ccf1)

    bis, int = GRASS.calc_bisector(v_grid, ccf1)
    bis_inv_slope = GRASS.calc_bisector_inverse_slope(bis, int)

    # plt.plot(mean(bis, dims=2), mean(int, dims=2), c=colors[i])

    plt.scatter(rvs1 .- mean(rvs1), bis_inv_slope .- mean(bis_inv_slope), c=colors[i], s=2)

end
# plt.xlabel(L"\Delta v\ {\rm (m/s)}")
# plt.ylabel(L"{\rm Flux}")
# plt.show()

plt.xlim(-1.5, 1.5)
plt.ylim(-1.0, 0.85)
plt.xlabel(L"\Delta {\rm RV\ (m/s)}")
plt.ylabel(L"{\rm BIS\ (m/s)}")
plt.show()
