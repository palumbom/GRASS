# environment + packages
using Pkg; Pkg.activate(".")
using CSV
using CUDA
using GRASS
using Peaks
using LsqFit
using Statistics
using DataFrames
using EchelleCCFs: λ_air_to_vac, calc_doppler_factor, MeasureRvFromCCFQuadratic as QuadraticFit

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

outdir = plotdir * "iag_comparison/"
if !isdir(outdir)
    mkdir(outdir)
end

# decide whether to use gpu
use_gpu = CUDA.functional()

# model the iag blends
function model_iag_blends(wavs_sim::AbstractArray{T,1}, flux_sim::AbstractArray{T,1},
                          wavs_iag::AbstractArray{T,1}, flux_iag::AbstractArray{T,1};
                          plot=false) where T<:Float64
    # calculate the residuals
    resids = flux_iag .- flux_sim

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

    # identify peaks in the residuals
    max_resid = 0.01
    minprom = 0.0005

    # guess the number of fits needed from minima
    m_inds = Peaks.argmaxima(.-resids, 10, strict=false)
    m_inds, m_proms = Peaks.peakproms(m_inds, .-resids, minprom=minprom, strict=false)
    m_inds, m_widths, m_left, m_right = Peaks.peakwidths(m_inds, .-resids, m_proms, strict=false)

    # convert width from pixels to wavelength
    m_widths .*= (wavs_iag[2] - wavs_iag[1])

    # set initial guess parameters
    nfits = length(m_inds)
    pgrid = zeros(nfits, 3)
    p0 = Array{Float64,1}[]
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

    return flux_iag./(tel_model(wavs_iag, fit.param)./flux_sim)
end

# figure 3 -- compare synthetic and IAG spectra + bisectors
function main()
    # get data
    lp = GRASS.LineProperties()
    files = GRASS.get_file(lp)

    # wavelength of line to synthesize/compare to iag
    for (i, file) in enumerate(files)
        # if !contains(file, "FeI_5434")
        #     continue
        # end

        # get the line name
        line_name = split(splitdir(file)[end], ".h5")[1]

        # get properties from line
        airwav = lp.λrest[i]
        depth = lp.depth[i]

        # get IAG spectrum and normalize it
        wavs_iag, flux_iag = GRASS.read_iag_atlas(isolate=true, airwav=airwav, buffer=1.5)
        flux_iag ./= maximum(flux_iag)

        # convolve IAG spectrum to LARS resolution
        wavs_iag, flux_iag = GRASS.convolve_gauss(wavs_iag, flux_iag, new_res=7e5, oversampling=5.0)

        # get depth from IAG spectrum
        idx1 = findfirst(x -> x .<= airwav - 0.125, wavs_iag)
        idx2 = findfirst(x -> x .>= airwav + 0.125, wavs_iag)
        botind = argmin(view(flux_iag, idx1:idx2)) + idx1
        iag_depth = 1.0 - minimum(view(flux_iag, idx1:idx2))

        # get velocity of IAG spectrum
        vels_iag, ccf_iag = calc_ccf(wavs_iag, flux_iag, [airwav], [iag_depth], mean((wavs_iag[2:end])./diff(wavs_iag)))
        rvs_iag, sigs_iag = calc_rvs_from_ccf(vels_iag, ccf_iag, fit_type=QuadraticFit)

        # shift the IAG spectrum to rest frame
        wavs_iag ./= calc_doppler_factor(rvs_iag)

        # set up for GRASS spectrum simulation
        # TODO fix depth issue!!
        function calculate_ideal_depth()
            dep_iter = true
            iters = 0
            while dep_iter
                # simulate the spectrum
                lines = [airwav]
                depths = [iag_depth - 0.05 + iters * 0.0025]
                templates = [file]
                resolution = 1e6
                spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution)
                disk = DiskParams(N=132, Nt=5)

                # simulate the spectrum
                wavs_sim, flux_sim = synthesize_spectra(spec, disk, use_gpu=use_gpu, verbose=false)
                flux_sim = dropdims(mean(flux_sim, dims=2), dims=2)

                # get the depth difference
                dep_diff = (1.0 - minimum(flux_sim)) - iag_depth
                if (abs(dep_diff) < 0.01) | (dep_diff > 0.0)
                    dep_iter = false
                    return depths[1]
                end
                iters += 1
            end
        end

        idepth = calculate_ideal_depth()

        # simulate the spectrum
        lines = [airwav]
        depths = [idepth]
        templates = [file]
        resolution = 1e6
        spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution)
        disk = DiskParams(N=132, Nt=5)

        # simulate the spectrum
        wavs_sim, flux_sim = synthesize_spectra(spec, disk, use_gpu=use_gpu)
        flux_sim = dropdims(mean(flux_sim, dims=2), dims=2)

        # get velocity of simulated spectrum
        vels_sim, ccf_sim = calc_ccf(wavs_sim, flux_sim, spec)
        rvs_sim, sigs_sim = calc_rvs_from_ccf(vels_sim, ccf_sim, fit_type=QuadraticFit)

        # shift the simulated spectrum to rest frame
        wavs_sim ./= calc_doppler_factor(rvs_sim)

        # interpolate IAG onto synthetic wavelength grid
        itp = GRASS.linear_interp(wavs_iag, flux_iag, bc=NaN)
        flux_iag = itp.(wavs_sim)
        wavs_iag = wavs_sim

        # clean the IAG spectrum
        flux_iag_cor = model_iag_blends(wavs_sim, flux_sim, wavs_iag, flux_iag)

        # get the synthetic + iag bisectors
        bis_sim, int_sim = GRASS.calc_bisector(wavs_sim, flux_sim)
        bis_iag, int_iag = GRASS.calc_bisector(wavs_iag, flux_iag)
        bis_iag_clean, int_iag_clean = GRASS.calc_bisector(wavs_iag, flux_iag_cor)

        # get ccfs
        v_grid_sim, ccf_sim = calc_ccf(wavs_sim, flux_sim, spec, normalize=true)
        v_grid_iag, ccf_iag = calc_ccf(wavs_iag, flux_iag, [wavs_iag[argmin(flux_iag)]],
                                       [1.0 - minimum(flux_iag)], 7e5, normalize=true)
        v_grid_iag_clean, ccf_iag_clean = calc_ccf(wavs_iag, flux_iag_cor, [wavs_iag[argmin(flux_iag_cor)]],
                                                   [1.0 - minimum(flux_iag_cor)], 7e5, normalize=true)

        # get ccf bisectors
        vel_sim, ccf_int_sim = GRASS.calc_bisector(v_grid_sim, ccf_sim)
        vel_iag, ccf_int_iag = GRASS.calc_bisector(v_grid_iag, ccf_iag)
        vel_iag_clean, ccf_int_iag_clean = GRASS.calc_bisector(v_grid_iag_clean, ccf_iag_clean)

        # align to arbitrary velocity
        vel_sim .-= (vel_sim[2])
        vel_iag .-= (vel_iag[2])
        vel_iag_clean .-= (vel_iag_clean[2])

        # big function for plotting
        function comparison_plots()
            # overplot the spectra
            fig = plt.figure()
            gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
            ax1 = fig.add_subplot(gs[1])
            ax2 = fig.add_subplot(gs[2])
            ax1.plot(wavs_sim, flux_sim, c="black", lw= 1.5, label=L"{\rm Synthetic}")
            ax1.plot(wavs_iag, flux_iag./maximum(flux_iag), marker="s", c="tab:blue", ms=2.0, lw=1.0, markevery=2, label=L"{\rm IAG}")
            ax1.plot(wavs_iag, flux_iag_cor./maximum(flux_iag_cor), alpha=0.9, marker="o", c="tab:green", ms=2.0, lw=1.0, markevery=2, label=L"{\rm Cleaned\ IAG}")
            ax2.plot(wavs_sim, flux_iag./maximum(flux_iag) .- flux_sim, c="tab:blue", marker="s", ms=2.0, lw=0)
            ax2.plot(wavs_sim, flux_iag_cor./maximum(flux_iag_cor) .- flux_sim, c="tab:green", marker="o", ms=2.0, lw=0)

            # set limits
            min_idx = argmin(flux_iag)
            ax1.set_xlim(airwav - 0.375, airwav + 0.375)
            ax1.set_ylim(minimum(flux_sim) - 0.1, 1.1)
            ax2.set_xlim(airwav - 0.5, airwav + 0.5)
            ax2.set_ylim(-0.125, 0.125)
            # ax2.set_ylim(minimum(flux_iag./maximum(flux_iag) .- flux_sim) - 0.04, maximum(flux_iag./maximum(flux_iag) .- flux_sim) + 0.04)

            # set tick labels, axis labels, etc.
            ax1.set_xticklabels([])
            ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
            ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
            ax2.set_ylabel(L"{\rm IAG\ -\ Synthetic}")
            ax1.legend()
            fig.tight_layout()

            # set the title
            ax1.set_title(("\${\\rm " * replace(line_name, "_" => "\\ ") * "}\$"))

            # save the plot
            fig.savefig(outdir * line_name * "_line.pdf")
            plt.clf(); plt.close()

            # plot the bisectors
            fig = plt.figure()
            gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
            ax1 = fig.add_subplot(gs[1])
            ax2 = fig.add_subplot(gs[2])

            ax1.plot(vel_sim[2:end-2], ccf_int_sim[2:end-2], color="black", lw=2.0, label=L"{\rm Synthetic}")
            ax1.plot(vel_iag[2:end-2], ccf_int_iag[2:end-2], marker="s", c="tab:blue", ms=2.0, lw=1.0, label=L"{\rm IAG}")
            ax1.plot(vel_iag_clean[2:end-2], ccf_int_iag_clean[2:end-2], marker="o", c="tab:green", ms=2.0, lw=1.0, label=L"{\rm Cleaned\ IAG}")
            ax2.plot(vel_iag[2:end-2] .- vel_sim[2:end-2], ccf_int_iag[2:end-2], c="tab:blue", marker="s", ms=2.0, lw=0.0)
            ax2.plot(vel_iag_clean[2:end-2] .- vel_sim[2:end-2], ccf_int_iag[2:end-2], c="tab:green", marker="o", ms=2.0, lw=0.0)

            # set tick labels, axis labels, etc.
            ax2.set_yticklabels([])
            ax2.yaxis.tick_right()
            # ax1.set_xlim(5434.4, 5434.6)
            ax1.set_ylim(0.1, 1.1)
            ax2.set_xlim(-20, 20)
            ax2.set_ylim(0.1, 1.1)
            ax1.set_xlabel(L"{\rm Relative\ Velocity\ (ms^{-1})}")
            ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
            ax2.set_xlabel(L"{\rm IAG\ -\ Synthetic\ (ms^{-1})}")
            ax1.legend(loc="upper right", prop=Dict("size"=>10), labelspacing=0.25)

            # set the title
            ax1.set_title(("\${\\rm " * replace(line_name, "_" => "\\ ") * "}\$"))

            # save the plot
            fig.savefig(outdir * line_name * "_bisector.pdf")
            plt.clf(); plt.close()
            return nothing
        end
        comparison_plots()
    end
    return nothing
end

if (run | plot)
    main()
end
