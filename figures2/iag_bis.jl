# environment + packages
using Pkg; Pkg.activate(".")
using CSV
using CUDA
using GRASS
using LsqFit
using Statistics
using DataFrames
using EchelleCCFs: λ_air_to_vac

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# decide whether to use gpu
use_gpu = CUDA.functional()

function interpolate_spec(wavs, flux)
    wavs_itp = collect(range(wavs[1], wavs[end], step=mean(diff(wavs))))
    itp = LinearInterpolation(wavs, flux)
    flux_itp = itp.(wavs_itp)
    return wavs_itp, flux_itp
end

# model the iag blends
function model_iag_blends(wavs_sim, flux_sim, wavs_iag, flux_iag; plot=false)
    # roughly align the spectra
    off1 = wavs_iag[argmin(flux_iag)] - wavs_sim[argmin(flux_sim)]
    wavs_iag .-= off1

    # interpolate onto the same grid
    flux_iag ./= maximum(flux_iag)
    itp = LinearInterpolation(wavs_sim, flux_sim, extrapolation_bc=1.0)
    wavs_sim = wavs_iag
    flux_sim = itp.(wavs_sim)

    # models for fit
    @. gaussian(x, a, b, c) = a * exp(-(x - b)^2/(2 * c^2)) + 1
    function tel_model(x, p)
        n = length(p) ÷ 3
        out = ones(length(x))
        for i in 1:n
            out .*= gaussian(x, p[3i-2:3i]...)
        end
        return out .* flux_sim
    end

    # calculate some residuals
    p0 = Array{Float64,1}[]
    resids = flux_iag .- flux_sim
    max_resid = 0.01
    minprom = 0.0025
    iters = 0
    while maximum(abs.(resids)) >= max_resid
        # iterate counter
        iters += 1

        # guess the number of fits needed from minima
        m_inds = Peaks.argmaxima(.-resids, 12, strict=false)
        m_inds, m_proms = Peaks.peakproms(m_inds, .-resids, minprom=minprom, strict=false)
        m_inds, m_widths, m_left, m_right = Peaks.peakwidths(m_inds, .-resids, m_proms, strict=false)

        # convert width from pixels to wavelength
        m_widths .*= (wavs_iag[2] - wavs_iag[1])

        # set initial guess parameters
        nfits = length(m_inds)
        pgrid = zeros(nfits, 3)
        for i in 1:nfits
            thresh = abs(wavs_sim[argmin(flux_sim)] - wavs_iag[m_inds[i]])
            if thresh <= 0.1
                continue
            end
            pgrid[i,1] = -m_proms[i]
            pgrid[i,2] = wavs_iag[m_inds[i]]
            pgrid[i,3] = m_widths[i]
        end
        p0 = Array{Float64,1}(vcat(p0, pgrid'...))

        # do the fit
        fit = curve_fit(tel_model, wavs_iag, flux_iag, p0)
        resids .= flux_iag./(tel_model(wavs_iag, fit.param)./flux_sim) .- flux_sim

        # plot diagnostics
        if plot
            fig = plt.figure(figsize=(8,6))
            gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
            ax1 = fig.add_subplot(gs[1])
            ax2 = fig.add_subplot(gs[2])
            ax1.plot(wavs_iag, flux_iag, label="Observed IAG", color="tab:blue")
            ax1.plot(wavs_iag, tel_model(wavs_iag, fit.param), label="Modeled IAG ", ls=":", color="tab:green")
            for i in m_inds
                thresh = abs(wavs_sim[argmin(flux_sim)] - wavs_iag[i])
                if thresh <= 0.1
                    continue
                end
                ax1.axvline(wavs_sim[i])
            end
            ax2.scatter(wavs_sim, flux_iag./tel_model(wavs_iag, fit.param), c="k", s=0.5)
            ax1.legend()
            ax1.set_xticklabels([])
            ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
            ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
            ax2.set_ylabel(L"{\rm IAG/Model}")
            plt.show()
            plt.clf(); plt.close()
        end

        # return condition
        if (iters > 5) || (maximum(abs.(resids)) < max_resid)
            return flux_iag./(tel_model(wavs_iag, fit.param)./flux_sim)
        end
    end
    return flux_iag
end

# figure 3 -- compare synthetic and IAG spectra + bisectors
function main()
    # get data
    lp = GRASS.LineProperties()
    files = GRASS.get_file(lp)

    # wavelength of line to synthesize/compare to iag
    for (i, file) in enumerate(files)
        if !contains(file, "FeI_5434")
            continue
        end

        # get properties from line
        airwav = lp.λrest[i]
        depth = lp.depth[i]

        # get IAG spectrum and normalize it
        wavs_iag, flux_iag = read_iag(isolate=true, airwav=airwav, buffer=2.0)
        flux_iag ./= maximum(flux_iag)

        # get depth from IAG spectrum
        idx1 = findfirst(x -> x .<= airwav - 0.25, wavs_iag)
        idx2 = findfirst(x -> x .>= airwav + 0.25, wavs_iag)
        botind = argmin(view(flux_iag, idx1:idx2)) + idx1
        iag_depth = 1.0 - minimum(view(flux_iag, idx1:idx2))

        # set up for GRASS spectrum simulation
        # TODO fix depth issue!!
        lines = [airwav]
        depths = [iag_depth + 0.0125]
        templates = [file]
        resolution = 1e6
        spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution)
        disk = DiskParams(N=132, Nt=10)

        # simulate the spectrum
        wavs_sim, flux_sim = synthesize_spectra(spec, disk, use_gpu=use_gpu)
        flux_sim = dropdims(mean(flux_sim, dims=2), dims=2)

        # get the synthetic bisector
        bis_sim, int_sim = GRASS.calc_bisector(wavs_sim, flux_sim)

        # align the spectra
        offset = wavs_sim[argmin(flux_sim)] - wavs_iag[botind]
        wavs_iag .+= offset

        # interpolate the IAG spectrum onto same grid as synthetic spectrum
        itp = GRASS.linear_interp(wavs_iag, flux_iag)
        wavs_iag_new = range(minimum(wavs_iag), maximum(wavs_iag), step=minimum(diff(wavs_sim)))
        flux_iag = itp.(wavs_iag_new)
        wavs_iag = wavs_iag_new

        # get the bisector
        bis_iag, int_iag = GRASS.calc_bisector(wavs_iag, flux_iag)

        println(minimum(flux_sim))
        println(minimum(flux_iag))


        plt.plot(wavs_iag, flux_iag)
        plt.plot(wavs_sim, flux_sim)
        plt.show()

        plt.plot(bis_iag, int_iag)
        plt.plot(bis_sim, int_sim)
        plt.show()

    end

    # airwav = 5434.5232;
    # geff = [0.0];

    # # get spectrum and interpolate onto even wavelength grid
    # wavs_iag, flux_iag = read_iag(airwav=airwav)
    # wavs_iag, flux_iag = interpolate_spec(wavs_iag, flux_iag)

    # # set up for GRASS spectrum simulation
    # lines = [airwav]
    # depths = [1.0 - minimum(flux_iag) + 0.05]
    # res = 7e5
    # spec = SpecParams(lines=lines, depths=depths, geffs=geff, resolution=res)
    # disk = DiskParams(N=132, Nt=15)

    # # synthesize spectra,
    # lambdas1, outspec1 = synthesize_spectra(spec, disk, seed_rng=false, verbose=true, top=NaN, use_gpu=use_gpu)

    # # calculate ccf, and get CCF bisector
    # v_grid, ccf1 = calc_ccf(lambdas1, outspec1, spec, normalize=true)
    # outspec1 = mean(outspec1, dims=2)[:,1]
    # ccfm = mean(ccf1, dims=2)[:,1]
    # vel_sim, bis_sim = GRASS.measure_bisector(v_grid, ccfm, interpolate=false, top=btop, len=len)

    # # model the line blends out of the IAG spectrum
    # println(">>> Modeling out line blends in IAG spectrum...")
    # clean = true
    # if clean
    #     flux_iag_cor = model_iag_blends(lambdas1, outspec1, wavs_iag, flux_iag)
    # else
    #     flux_iag_cor = flux_iag
    # end

    # # get offsets to align the spectra in wavelength
    # off1 = wavs_iag[argmin(flux_iag)] - lambdas1[argmin(outspec1[:,1])]
    # off2 = wavs_iag[argmin(flux_iag_cor)+1] - lambdas1[argmin(outspec1[:,1])]
    # wavs_iag .-= off1

    # # calculate a CCF for the IAG spectrum and trim it
    # v_grid_iag, ccf_iag = calc_ccf(wavs_iag, flux_iag, [wavs_iag[argmin(flux_iag)]],
    #                                [1.0 - minimum(flux_iag)],
    #                                1e6, normalize=true)

    # ind1 = findfirst(x -> x .> -vlim, v_grid_iag)
    # ind2 = length(v_grid_iag)#findfirst(x -> x .> vlim, v_grid_iag)
    # vel_iag, bis_iag = GRASS.measure_bisector(v_grid_iag[ind1:ind2], ccf_iag[ind1:ind2],
    #                                           interpolate=true, top=btop, len=len)

    # # calculate a CCF for the cleaned IAG spectrum and trim it
    # wavs_iag, flux_iag_cor = interpolate_spec(wavs_iag, flux_iag_cor)
    # v_grid_iag2, ccf_iag2 = calc_ccf(wavs_iag, flux_iag_cor, [wavs_iag[argmin(flux_iag_cor)]],
    #                                 [1.0 - minimum(flux_iag_cor)],
    #                                 1e6, normalize=true)

    # ind1 = findfirst(x -> x .> -vlim, v_grid_iag2)
    # ind2 = findfirst(x -> x .> vlim, v_grid_iag2)
    # ind2 = length(v_grid_iag2)
    # vel_iag2, bis_iag2 = GRASS.measure_bisector(v_grid_iag2[ind1:ind2], ccf_iag2[ind1:ind2],
    #                                             interpolate=true, top=btop, len=len)

    # # interpolate IAG onto same wavelength scale as synthetic spectrum
    # itp = LinearInterpolation(wavs_iag, flux_iag, extrapolation_bc=1.0)
    # flux_iag_itp = itp.(lambdas1)
    # itp = LinearInterpolation(wavs_iag, flux_iag_cor, extrapolation_bc=1.0)
    # flux_iag_itp2 = itp.(lambdas1)

    # function comparison_plots()
    #     # overplot the spectra
    #     fig = plt.figure()
    #     gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
    #     ax1 = fig.add_subplot(gs[1])
    #     ax2 = fig.add_subplot(gs[2])
    #     ax1.plot(lambdas1, outspec1, c="black", lw= 1.5, label=L"{\rm Synthetic}")
    #     ax1.plot(wavs_iag, flux_iag./maximum(flux_iag), marker="s", c="tab:blue", ms=2.0, lw=1.0, markevery=10, label=L"{\rm IAG}")
    #     ax1.plot(wavs_iag, flux_iag_cor./maximum(flux_iag_cor), alpha=0.9, marker="o", c="tab:green", ms=2.0, lw=1.0, markevery=10, label=L"{\rm Cleaned\ IAG}")
    #     ax2.plot(lambdas1, flux_iag_itp./maximum(flux_iag_itp) .- outspec1, c="tab:blue", marker="s", ms=2.0, lw=0)
    #     ax2.plot(lambdas1, flux_iag_itp2./maximum(flux_iag_itp2) .- outspec1, c="tab:green", marker="o", ms=2.0, lw=0)

    #     # set tick labels, axis labels, etc.
    #     ax1.set_xticklabels([])
    #     ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
    #     ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
    #     ax2.set_ylabel(L"{\rm IAG\ -\ Synthetic}")
    #     ax1.legend()
    #     fig.tight_layout()

    #     # save the plot
    #     fig.savefig(plotdir * "fig3a.pdf")
    #     plt.clf(); plt.close()
    #     println(">>> Figure written to: " * plotdir * "fig3a.pdf")

    #     # align bisectors to arbitrary point
    #     vel_sim .-= mean(vel_sim) #.- 42.0
    #     vel_iag .-= mean(vel_iag) #.- 42.0
    #     vel_iag2 .-=  mean(vel_iag2) #.- 42.0

    #     # plot the bisectors
    #     fig = plt.figure()
    #     gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
    #     ax1 = fig.add_subplot(gs[1])
    #     ax2 = fig.add_subplot(gs[2])
    #     ax1.plot(vel_sim[2:end-2], bis_sim[2:end-2], color="black", lw=2.0, label=L"{\rm Synthetic}")
    #     ax1.plot(vel_iag[2:end-2], bis_iag[2:end-2], marker="s", c="tab:blue", ms=2.0, lw=1.0, label=L"{\rm IAG}")
    #     ax1.plot(vel_iag2[2:end-2], bis_iag2[2:end-2], marker="o", c="tab:green", ms=2.0, lw=1.0, label=L"{\rm Cleaned\ IAG}")
    #     ax2.plot(vel_iag[2:end-2] .- vel_sim[2:end-2], bis_iag[2:end-2], c="tab:blue", marker="s", ms=2.0, lw=0.0)
    #     ax2.plot(vel_iag2[2:end-2] .- vel_sim[2:end-2], bis_iag[2:end-2], c="tab:green", marker="o", ms=2.0, lw=0.0)

    #     # set tick labels, axis labels, etc.
    #     ax2.set_yticklabels([])
    #     ax2.yaxis.tick_right()
    #     ax1.set_ylim(0.1, 1.1)
    #     ax2.set_xlim(-20, 20)
    #     ax2.set_ylim(0.1, 1.1)
    #     ax1.set_xlabel(L"{\rm Relative\ Velocity\ (ms^{-1})}")
    #     ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
    #     ax2.set_xlabel(L"{\rm IAG\ -\ Synthetic\ (ms^{-1})}")
    #     ax1.legend(loc="upper right", prop=Dict("size"=>10), labelspacing=0.25)

    #     # save the plot
    #     fig.savefig(plotdir * "fig3b.pdf")
    #     plt.clf(); plt.close()
    #     println(">>> Figure written to: " * plotdir * "fig3b.pdf")
    #     return nothing
    # end

    # comparison_plots()
    return nothing

    return nothing
end

if (run | plot)
    main()
end
