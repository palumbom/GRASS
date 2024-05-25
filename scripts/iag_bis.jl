# environment + packages
using CSV
using CUDA
using GRASS
using Peaks
using Optim
using LsqFit
using SPICE
using Statistics
using DataFrames
using EchelleCCFs
using EchelleCCFs: λ_air_to_vac, calc_doppler_factor, MeasureRvFromCCFQuadratic as QuadraticFit

GRASS.get_kernels()

neid_timestamps = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
#convert from utc to et as needed by SPICE
time_stamps = utc2et.(neid_timestamps)
#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

# plotting
# using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
# mpl.style.use(GRASS.moddir * "fig.mplstyle")
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

# get command line args and output directories
# include(joinpath(abspath(@__DIR__), "paths.jl"))
plotdir = string(abspath(joinpath("figures", "iag_comparison_NEIDeclipse")))
datadir = string(abspath("data"))

if !isdir(plotdir)
    mkdir(plotdir)
end

# decide whether to use gpu
use_gpu = false
# @assert CUDA.functional()

# get data
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
files = GRASS.get_file(lp)
line_names = GRASS.get_name(lp)

# read in optimized depths
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

df_tuned = CSV.read(joinpath(datadir, "tuned_params.csv"), DataFrame)

# wavelength of line to synthesize/compare to iag
for (i, file) in enumerate(files)
    # if !contains(file, "FeI_5383")
    #     continue
    # end
    println(">>> Running " * line_names[i] * "...")

    # i = 5
    # file = line_names[i]

    # get properties from line
    line_name = line_names[i]
    airwav = lp.λrest[i]
    depth = lp.depth[i]

    # get IAG spectrum and normalize it
    wavs_iag0, flux_iag0 = GRASS.read_iag_atlas(isolate=true, airwav=airwav, buffer=1.5)
    flux_iag0 ./= maximum(flux_iag0)

    # convolve IAG spectrum to LARS resolution
    wavs_iag, flux_iag = GRASS.convolve_gauss(wavs_iag0, flux_iag0, new_res=1e6, oversampling=4.0)

   # get depth from IAG spectrum
    buff = 0.12575
    if contains("FeI_5383", line_name)
        buff = 0.3
    elseif contains("FeI_5434", line_name)
        buff = 0.3
    elseif contains("FeI_5382", line_name)
        buff = 0.2
    elseif contains("FeI_5576", line_name)
        buff = 0.25
    elseif contains("CaI_6169.0", line_name)
        buff = 0.25
    elseif contains("CaI_6169.5", line_name)
        buff = 0.1475
    elseif contains("FeI_6170", line_name)
        buff = 0.175
    elseif contains("FeI_6301", line_name)
        buff = 0.25
    elseif contains("FeI_6302", line_name)
        buff = 0.125
    end

    idxl = findfirst(x -> x .>= airwav - buff, wavs_iag)
    idxr = findfirst(x -> x .>= airwav + buff, wavs_iag)
    iag_bot = minimum(view(flux_iag, idxl:idxr))
    iag_depth = 1.0 - iag_bot

    # get the depth for the simulation
    sim_depth = df[i, "optimized_depth"]

    # simulate the spectrum
    lines = [airwav]
    depths = [sim_depth]
    templates = [file]
    resolution = 7e5
    spec = SpecParams(lines=lines, depths=depths, templates=templates,
                      resolution=resolution, buffer=1.5, oversampling=2.0)
    N = 50
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    # simulate the spectrum
    wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, obs_long, obs_lat, alt, lines ./ 10.0, time_stamps, verbose=true, use_gpu=false)
    flux_sim = dropdims(mean(flux_sim, dims=2), dims=2)

    # interpolate iag on synth wavelength grid
    itp = GRASS.linear_interp(wavs_iag, flux_iag, bc=NaN)
    flux_iag = itp.(wavs_sim)
    wavs_iag = copy(wavs_sim)

    # get width in velocity for CCF
    idxl_sim, idxr_sim = GRASS.find_wing_index(0.95, flux_sim)
    if contains("FeI_6302", line_name)
        idxl_sim, idxr_sim = GRASS.find_wing_index(0.875, flux_sim)
    end

    # get width in angstroms
    width_ang = wavs_sim[idxr_sim] - wavs_sim[idxl_sim]

    # convert to velocity
    width_vel = GRASS.c_ms * width_ang / wavs_sim[argmin(flux_sim)]
    Δv_max = round((width_vel + 1e3)/100) * 100

    # calculate ccfs for spectra
    if contains("FeI_6302", line_name)
        Δv_max -= 200
    end
    v_grid_iag, ccf_iag = GRASS.calc_ccf(wavs_iag, flux_iag, lines, depths,
                                         7e5, Δv_step=100.0, Δv_max=Δv_max,
                                         mask_type=EchelleCCFs.GaussianCCFMask)

    v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths,
                                         7e5, Δv_step=100.0, Δv_max=Δv_max,
                                         mask_type=EchelleCCFs.GaussianCCFMask)

    # deal with annoying line blend
    if contains("FeI_6302", line_name)
        amin_ccf = argmin(ccf_iag)
        idx_b = findfirst(x -> x .> 0.925, ccf_iag[amin_ccf:end]) + amin_ccf
        ccf_iag = ccf_iag[1:idx_b]
        v_grid_iag = v_grid_iag[1:idx_b]
    end

    # plt.plot(wavs_iag, flux_iag)
    # plt.plot(wavs_sim, flux_sim)
    # plt.axvline(v_grid_iag[idx_b])
    # plt.plot(v_grid_iag, ccf_iag)
    # plt.plot(v_grid_sim, ccf_sim)
    # plt.show()
    # break

    # get bisectors
    top = 0.9
    vel_iag, int_iag = GRASS.calc_bisector(v_grid_iag, ccf_iag, nflux=50, top=top)
    vel_sim, int_sim = GRASS.calc_bisector(v_grid_sim, ccf_sim, nflux=50, top=top)

    # plt.plot(vel_iag, int_iag)
    # plt.plot(vel_sim, int_sim)
    # plt.show()

    # compute velocity as mean bisector between N and M % depth
    N = 0.20
    M = 0.70
    idx1 = findfirst(x -> x .>= N * iag_depth + iag_bot, int_iag)
    idx2 = findfirst(x -> x .>= M * iag_depth + iag_bot, int_iag)
    if isnothing(idx2)
        idx2 = findfirst(x -> x .>= 0.9, int_iag)
    end
    rv_iag = mean(view(vel_iag, idx1:idx2))

    idx1 = findfirst(x -> x .>= N * sim_depth + minimum(flux_sim), int_sim)
    idx2 = findfirst(x -> x .>= M * sim_depth + minimum(flux_sim), int_sim)
    if isnothing(idx2)
        idx2 = findfirst(x -> x .>= 0.9, int_sim)
    end
    rv_sim = mean(view(vel_sim, idx1:idx2))

    # transform to lab frame
    vel_iag .-= rv_iag
    vel_sim .-= rv_sim
    wavs_iag ./= calc_doppler_factor(rv_iag)
    wavs_sim ./= calc_doppler_factor(rv_sim)

    # interpolate IAG onto synthetic wavelength grid
    itp = GRASS.linear_interp(wavs_iag, flux_iag)
    flux_iag = itp.(wavs_sim)
    wavs_iag = copy(wavs_sim)

    # re-compute line isolation indices because of interpolation
    idxl = findfirst(x -> x .>= airwav - buff, wavs_iag)
    idxr = findfirst(x -> x .>= airwav + buff, wavs_iag)

    # recompute bisectors b/c of interpolation
    v_grid_iag, ccf_iag = GRASS.calc_ccf(wavs_iag, flux_iag, lines, depths,
                                         7e5, Δv_step=100.0, Δv_max=Δv_max,
                                         mask_type=EchelleCCFs.GaussianCCFMask)

    v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths,
                                         7e5, Δv_step=100.0, Δv_max=Δv_max,
                                         mask_type=EchelleCCFs.GaussianCCFMask)

    if contains("FeI_6302", line_name)
        amin_ccf = argmin(ccf_iag)
        idx_b = findfirst(x -> x .> 0.925, ccf_iag[amin_ccf:end]) + amin_ccf
        ccf_iag = ccf_iag[1:idx_b]
        v_grid_iag = v_grid_iag[1:idx_b]
    end

    # # set errant ccf values
    # idx0 = iszero.(ccf_iag)
    # idxnz = findfirst(x -> x .> 0.0, ccf_iag)
    # ccf_iag[idx0] .= ccf_iag[idxnz]

    # get bisectors
    vel_iag, int_iag = GRASS.calc_bisector(v_grid_iag, ccf_iag, nflux=50, top=top)
    vel_sim, int_sim = GRASS.calc_bisector(v_grid_sim, ccf_sim, nflux=50, top=top)

    # find mean velocities in order to align bisectors
    N = 0.20
    M = 0.70
    idx1 = findfirst(x -> x .>= N * sim_depth + minimum(flux_sim), int_sim)
    idx2 = findfirst(x -> x .>= M * sim_depth + minimum(flux_sim), int_sim)
    if isnothing(idx2)
        idx2 = findfirst(x -> x .>= 0.9, int_sim)
    end
    rv_sim = mean(view(vel_sim, idx1:idx2))

    idx1 = findfirst(x -> x .>= N * iag_depth + iag_bot, int_iag)
    idx2 = findfirst(x -> x .>= M * iag_depth + iag_bot, int_iag)
    if isnothing(idx2)
        idx2 = findfirst(x -> x .>= 0.9, int_iag)
    end
    rv_iag = mean(view(vel_iag, idx1:idx2))

    # align the bisectors
    vel_sim .-= rv_sim
    vel_iag .-= rv_iag

    # get the bisector residuals
    bis_resids = vel_iag .- vel_sim
    med_resids = median(bis_resids)

    vel_sim .+= med_resids

    # get the tuned fluxes for BIS
    b1 = df_tuned[i, "b1"]
    b2 = df_tuned[i, "b2"]
    b3 = df_tuned[i, "b3"]
    b4 = df_tuned[i, "b4"]

    # convert to flux level
    dep_sim = 1.0 - minimum(ccf_sim)
    i1 = 1.0 - b1 * dep_sim
    i2 = 1.0 - b2 * dep_sim
    i3 = 1.0 - b3 * dep_sim
    i4 = 1.0 - b4 * dep_sim

    # big function for plotting
    function comparison_plots()
        # make plot objects
        fig = plt.figure(figsize=(6.4,4.8))
        gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
        ax1 = fig.add_subplot(gs[1])
        ax2 = fig.add_subplot(gs[2])

        # plot spectra
        ax1.plot(wavs_sim, flux_sim, marker="o", c="black", ms=3.0, lw=1.5, markevery=2, label="Synthetic")
        ax1.plot(wavs_iag, flux_iag, marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=2, label="IAG")
        # ax1.plot(wavs_iag, flux_mod, alpha=0.9, marker="o", c=colors[2], ms=2.0, lw=1.0, markevery=2, label=L"{\rm Cleaned\ IAG}")

        # plot resids
        ax2.plot(wavs_sim, flux_iag .- flux_sim, c=colors[1], marker="s", ms=2.0, lw=0, markevery=2)
        # ax2.plot(wavs_sim, flux_mod .- flux_sim, c=colors[2], marker="^", ms=2.0, lw=0, markevery=2)

        # find limits
        idx_min = argmin(flux_sim)
        idx1 = idx_min - findfirst(x -> x .> 0.95, flux_sim[idx_min:-1:1])
        idx2 = idx_min + findfirst(x -> x .> 0.95, flux_sim[idx_min:end])

        # set limits
        min_idx = argmin(flux_iag)
        ax1.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
        ax1.set_ylim(minimum(flux_sim) - 0.1, 1.1)
        ax2.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
        ax2.set_ylim(-0.0575, 0.0575)

        # set tick labels, axis labels, etc.
        ax1.set_xticklabels([])
        # ax1.set_ylabel(L"{\rm Normalized\ Flux}", labelpad=15)
        # ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
        # ax2.set_ylabel(L"{\rm IAG\ -\ GRASS}")
        ax1.legend()
        # fig.tight_layout()

        # set the title
        # title = replace(line_name, "_" => "\\ ")
        # idx = findfirst('I', title)
        # title = title[1:idx-1] * "\\ " * title[idx:end] * "\\ \\AA"
        # ax1.set_title(("\${\\rm " * title * "}\$"))

        # save the plot
        fig.subplots_adjust(wspace=0.05)
        fig.savefig(joinpath(plotdir, line_name * "_line.pdf"))
        plt.clf(); plt.close()

        # plot the bisectors
        fig = plt.figure(figsize=(6.4,4.8))
        gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
        ax1 = fig.add_subplot(gs[1])
        ax2 = fig.add_subplot(gs[2])

        # plot bisectors
        ax1.plot(vel_sim[4:end], int_sim[4:end], marker="o", color="black", ms=3.0, lw=2.0, markevery=1, label="Synthetic")
        ax1.plot(vel_iag[4:end], int_iag[4:end], marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=1, label="IAG")
        # ax1.plot(vel_sim2, int_sim2, marker="o", color="black", ms=3.0, lw=2.0, markevery=1, label=L"{\rm Derp}")
        # ax1.plot(vel_mod[2:end], int_mod[2:end], marker="^", c=colors[2], ms=2.0, lw=1.0, markevery=1, label=L"{\rm Cleaned\ IAG}")

        # plot BIS levels
        # ax1.axhline(i1, ls="--", c="k")
        # ax1.axhline(i2, ls="--", c="k")
        # ax1.axhline(i3, ls=":", c="k")
        # ax1.axhline(i4, ls=":", c="k")

        # plot residuals
        ax2.plot(vel_iag[4:end] .- vel_sim[4:end], int_iag[4:end], c=colors[1], marker="s", ms=2.0, lw=0.0, markevery=1)
        # ax2.plot(vel_mod[2:end] .- vel_sim[2:end], int_mod[2:end], c=colors[2], marker="o", ms=2.0, lw=0.0, markevery=2)

        # set tick labels, axis labels, etc.
        ax2.set_yticklabels([])
        ax2.yaxis.tick_right()
        ax1.set_ylim(minimum(flux_sim) - 0.05, 1.05)
        ax2.set_xlim(-35, 35)
        ax2.set_ylim(minimum(flux_sim) - 0.05, 1.05)
        # ax1.set_xlabel(L"{\rm Relative\ Velocity\ (m\ s^{-1})}", fontsize=13)
        # ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
        # ax2.set_xlabel(L"{\rm IAG\ -\ GRASS\ (m\ s^{-1})}", fontsize=13)
        ax1.legend(labelspacing=0.25)
        # ax2.legend(loc="upper right", prop=Dict("size"=>12.5), labelspacing=0.25)

        # set the title
        # fig.suptitle(("\${\\rm " * title * "}\$"), y=0.95)

        # save the plot
        fig.subplots_adjust(hspace=0.05)
        fig.savefig(joinpath(plotdir, line_name * "_bisector.pdf"))
        plt.clf(); plt.close()
        return nothing
    end
    comparison_plots()
end
