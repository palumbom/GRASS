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
line_depths = GRASS.get_depth(lp)
line_names = GRASS.get_name(lp)
line_titles = replace.(line_names, "_" => " ")

for (idx, file) in enumerate(lp.file)
    # set up paramaters for spectrum
    N = 132
    Nt = 1250
    lines = [rest_wavelengths[idx]]
    depths = [line_depths[idx]]
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
    fig1, axs1 = plt.subplots(figsize=(10.0,7.7))
    fig2, axs2 = plt.subplots(figsize=(10.0,7.7), nrows=2, ncols=2, sharex=true, sharey=true)
    # fig3, axs3 = plt.subplots(figsize=(11,8.5))

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

        # set mask with to two pixels in velocity
        mask_width = (GRASS.c_ms/resolutions[i]) * (2.0/oversampling)

        # calculate a ccf
        v_grid, ccf1 = calc_ccf(wavs_out[:,1], flux_out, lines, depths,
                                resolutions[i], mask_width=mask_width,
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
        fit_label = "\$ " .* string(slope) .* "\$"

        # get the bisectors and subtract off mean
        mean_bis = dropdims(mean(bis, dims=2), dims=2)[2:end-1]
        mean_int = dropdims(mean(int, dims=2), dims=2)[2:end-1]
        mean_bis .-= mean(mean_bis)

        # plot the bisector
        axs1.plot(mean_bis, mean_int, c=colors[i], label=labels[i], lw=3.0)

        # plot BIS and apparent RV
        axs2[i].scatter(xdata, ydata, c=colors[i], s=2)#, label=labels[i])
        axs2[i].plot(xmodel, ymodel, c="k", ls="--", label=L"{\rm Slope } \approx\ " * fit_label, lw=2.5)
        axs2[i].text(-1.1, -1.05, labels[i], fontsize=18, bbox=Dict("boxstyle" => "round", "fc" => "white", "ec" => "lightgray"))

        # plot the line profiles
        # axs3.plot(mean(wavs_out, dims=2), mean(flux_out, dims=2), c=colors[i], label=labels[i])

        # plot ccfs
        # plt.plot(v_grid, mean(ccf1, dims=2), c=colors[i])
    end
    # set font sizes
    title_font = 26
    tick_font = 20
    legend_font = 18

    # set plot stuff for first plot
    axs1.set_xlabel(L"\Delta v\ {\rm (m s}^{-1}{\rm )}", fontsize=title_font)
    axs1.set_ylabel(L"{\rm CCF}", fontsize=title_font)
    axs1.set_title("\${\\rm " * replace(line_titles[idx], " "=>"\\ ") * "\\ \\AA}\$", fontsize=title_font)
    axs1.tick_params(axis="both", which="major", labelsize=tick_font+4)
    axs1.legend(fontsize=legend_font+2)
    fig1.savefig("plottos/" * line_names[idx] * "_bisector.pdf", bbox_inches="tight")

    # set plot stuff for second plot
    for ax in axs2
        ax.tick_params(axis="both", which="major", labelsize=tick_font)
        ax.legend(loc="upper right", fontsize=legend_font, ncols=1,
                  columnspacing=0.8, handletextpad=0.5, labelspacing=0.08)
        ax.set_xlim(-1.25, 1.25)
        ax.set_ylim(-1.25, 1.4)
        ax.set_xticks([-1, -0.5, 0.0, 0.5, 1.0])
        ax.set_yticks([-1, -0.5, 0.0, 0.5, 1.0])
    end
    fig2.supxlabel(L"\Delta v\ {\rm (m s}^{-1}{\rm )}", fontsize=title_font, x=0.55, y=0.05)
    fig2.supylabel(L"{\rm BIS}\ - \overline{\rm BIS}\ {\rm (m s}^{-1}{\rm )}", fontsize=title_font, x=0.03, y=0.52)
    fig2.suptitle("\${\\rm " * replace(line_titles[idx], " "=>"\\ ") * "\\ \\AA}\$", fontsize=title_font, x=0.55, y=0.95)
    plt.tight_layout(pad=-2.0)
    fig2.savefig("plottos/" * line_names[idx] * "_rv_vs_bis.pdf", bbox_inches="tight")

    # set plot stuff for third plot
    # axs3.tick_params(axis="both", which="major", labelsize=tick_font)
    # axs3.set_xlabel(L"{\rm Wavelength\ (\AA)}", fontsize=title_font)
    # axs3.set_ylabel(L"{\rm Normalized\ Flux}", fontsize=title_font)
    # axs3.set_title("\${\\rm " * replace(line_titles[idx], " "=>"\\ ") * "}\$", fontsize=title_font)
    # axs3.legend(fontsize=legend_font)
    # fig3.savefig("plottos/" * line_names[idx] * "_line_profile.pdf", bbox_inches="tight")

    plt.close("all")
end
