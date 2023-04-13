using Pkg; Pkg.activate(".")
using CUDA
using GRASS
using Printf
using Revise
using Statistics
using EchelleCCFs
using Polynomials
using BenchmarkTools
using HypothesisTests

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")

function round_and_format(num::Float64)
    rounded_num = Int(round(num))
    formatted_num = collect(string(rounded_num))
    num_length = length(formatted_num)

    if num_length <= 3
        return prod(formatted_num)
    end

    comma_idx = mod(num_length, 3)
    if comma_idx == 0
        comma_idx = 3
    end

    while comma_idx < num_length
        insert!(formatted_num, comma_idx+1, ',')
        comma_idx += 4
        num_length += 1
    end

    return replace(prod(formatted_num), "," => "{,}")
end

# get line properties
lp = GRASS.LineProperties()
line_species = GRASS.get_species(lp)
rest_wavelengths = GRASS.get_rest_wavelength(lp)
line_names = GRASS.get_name(lp)
line_titles = replace.(line_names, "_" => " ")

for (idx, file) in enumerate(lp.file)
    if line_names[idx] != "FeI_5576"
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
    lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=true,
                                            verbose=true, use_gpu=true)

    # set oversampling and resolutions, set color sequence
    oversampling = 4.0
    resolutions = reverse([0.98e5, 1.2e5, 1.37e5, 2.7e5])
    instruments = ["PEPSI", "EXPRES", "NEID", "KPF"]
    # colors = ["tab:blue", "tab:orange", "tab:green", "tab:purple"]
    colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

    # get labels for resolutions
    labels = "\$ R \\sim " .* round_and_format.(resolutions) .* "{\\rm \\ (" .* instruments .* ")}" .* "\$"

    # create fig + axes objects
    fig1, axs1 = plt.subplots(figsize=(11,8.5))
    fig2, axs2 = plt.subplots(figsize=(11,8.5), nrows=2, ncols=2, sharex=true, sharey=true)
    fig3, axs3 = plt.subplots(figsize=(11,8.5))

    # re-order axs so that its indexed by row
    axs2 = [axs2[1], axs2[3], axs2[2], axs2[4]]

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

        # calculate a ccf
        v_grid, ccf1 = calc_ccf(wavs_out[:,1], flux_out, lines, depths,
                                resolutions[i], normalize=true,
                                mask_type=EchelleCCFs.TopHatCCFMask,
                                Δv_step=125.0)
        rvs1, sigs1 = calc_rvs_from_ccf(v_grid, ccf1)

        # calculate bisector and BIS
        bis, int = GRASS.calc_bisector(v_grid, ccf1, nflux=20)
        bis_inv_slope = GRASS.calc_bisector_inverse_slope(bis, int)

        # subtract off mean
        xdata = rvs1 .- mean(rvs1)
        ydata = bis_inv_slope .- mean(bis_inv_slope)

        # fit to the BIS
        pfit = Polynomials.fit(xdata, ydata, 1)
        xmodel = range(-1.15, 1.15, length=5)
        ymodel = pfit.(xmodel)

        # get the slope of the fit
        slope = round(coeffs(pfit)[2], digits=3)
        fit_label = "\$ R \\sim " .* string(slope) .* "\$"

        # plot the bisectors
        axs1.plot(mean(bis, dims=2)[2:end-1], mean(int, dims=2)[2:end-1], c=colors[i], label=labels[i])

        # plot BIS and apparent RV
        axs2[i].scatter(xdata, ydata, c=colors[i], s=2, label=labels[i])
        axs2[i].plot(xmodel, ymodel, c="k", ls="--", label=L"{\rm Slope } \approx\ " * fit_label)

        # plot the line profiles
        axs3.plot(mean(wavs_out, dims=2), mean(flux_out, dims=2), c=colors[i], label=labels[i])

        # plot ccfs
        # plt.plot(v_grid, mean(ccf1, dims=2), c=colors[i])
    end
    # set plot stuff for first plot
    axs1.set_xlabel(L"\Delta v\ {\rm (m s}^{-1}{\rm )}")
    axs1.set_ylabel(L"{\rm CCF}")
    axs1.set_title("\${\\rm " * replace(line_titles[idx], " "=>"\\ ") * "}\$")
    axs1.legend()
    fig1.savefig("plottos/" * line_names[idx] * "_bisector.pdf", bbox_inches="tight")

    # set plot stuff for second plot
    for ax in axs2
        ax.legend(loc="upper right")
        ax.set_xlim(-1.25, 1.25)
        ax.set_ylim(-1.25, 1.25)
    end
    fig2.supxlabel(L"\Delta v\ {\rm (m s}^{-1}{\rm )}")
    fig2.supylabel(L"{\rm BIS}\ - \overline{\rm BIS}\ {\rm (m s}^{-1}{\rm )}")
    fig2.suptitle("\${\\rm " * replace(line_titles[idx], " "=>"\\ ") * "}\$")
    fig2.savefig("plottos/" * line_names[idx] * "_rv_vs_bis.pdf", bbox_inches="tight")

    # set plot stuff for third plot
    axs3.set_xlabel(L"{\rm Wavelength\ (\AA)}")
    axs3.set_ylabel(L"{\rm Normalized\ Flux}")
    axs3.set_title("\${\\rm " * replace(line_titles[idx], " "=>"\\ ") * "}\$")
    axs3.legend()
    fig3.savefig("plottos/" * line_names[idx] * "_line_profile.pdf", bbox_inches="tight")

    plt.close("all")
end
