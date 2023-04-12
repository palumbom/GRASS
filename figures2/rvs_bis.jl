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

# get line properties
lp = GRASS.LineProperties()
line_species = GRASS.get_species(lp)
rest_wavelengths = GRASS.get_rest_wavelength(lp)
line_names = GRASS.get_name(lp)
line_titles = replace.(line_names, "_" => " ")

for (idx, file) in enumerate(lp.file)
    if idx != 8
        continue
    end
    # set up paramaters for spectrum
    N = 132
    Nt = 1000
    lines = [rest_wavelengths[idx]]
    depths = [0.9]
    geffs = [0.0]
    templates = [file]
    variability = repeat([true], length(lines))
    resolution = 7e5
    seed_rng = true

    disk = DiskParams(N=N, Nt=Nt)
    spec1 = SpecParams(lines=lines, depths=depths, variability=variability,
                       geffs=geffs, templates=templates, resolution=resolution)
    lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=false, verbose=true, use_gpu=true)

    # set oversampling and resolutions, set color sequence
    oversampling = 2.0
    resolutions = reverse([0.98e5, 1.2e5, 1.37e5, 2.7e5])
    instruments = ["PEPSI", "EXPRES", "NEID", "KPF"]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:purple"]
    labels = "\$ R \\sim " .* string.(resolutions) .* "\$"

    # create fig + axes objects
    fig1, axs1 = plt.subplots()
    fig2, axs2 = plt.subplots(nrows=2, ncols=2, sharex=true, sharey=true)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    # loop over resolutions
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
                                mask_type=EchelleCCFs.TopHatCCFMask,
                                Δv_step=125.0)
        rvs1, sigs1 = calc_rvs_from_ccf(v_grid, ccf1)

        bis, int = GRASS.calc_bisector(v_grid, ccf1, nflux=20)
        bis_inv_slope = GRASS.calc_bisector_inverse_slope(bis, int)

        # plot the line profiles
        # axs1.plot(mean(wavs_out, dims=2), mean(flux_out, dims=2), c=colors[i], label=labels[i])

        # plot the bisectors
        # axs1.plot(mean(bis, dims=2)[2:end-1], mean(int, dims=2)[2:end-1], c=colors[i])

        # plot ccfs
        # plt.plot(v_grid, mean(ccf1, dims=2), c=colors[i])

        # plot BIS and apparent RV
        axs2[i].scatter(rvs1 .- mean(rvs1), bis_inv_slope .- mean(bis_inv_slope), c=colors[i], s=2)

    end
    # set plot stuff for first plot
    axs1.set_xlabel(L"\Delta v\ {\rm (m s}^{-1}{\rm )}")
    axs1.set_ylabel(L"{\rm Normalized\ Flux}")
    axs1.set_title(("\${\\rm " * replace(line_titles[idx], " "=>"\\ ") * "}\$"))
    axs1.legend()

    plt.show()

    # plt.xlim(-1.5, 1.5)
    # plt.ylim(-1.0, 0.85)
    # plt.xlabel(L"\Delta {\rm RV\ (m/s)}")
    # plt.ylabel(L"{\rm BIS\ (m/s)}")
    # plt.show()
end
