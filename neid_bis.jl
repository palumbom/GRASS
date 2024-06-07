# environment + packages
using CSV
using CUDA
using GRASS
using Peaks
using Optim
using LsqFit
using SPICE
using FITSIO
using FileIO
using JLD2
using Statistics
using DataFrames
using EchelleCCFs
using EchelleCCFs: λ_air_to_vac, calc_doppler_factor, MeasureRvFromCCFQuadratic as QuadraticFit
using RvSpectMLBase
using EchelleInstruments

GRASS.get_kernels()

neid_timestamps = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
#convert from utc to et as needed by SPICE
time_stamps = utc2et.(neid_timestamps)
#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

refractive_index = 1.0003

# plotting
# using LaTeXStrings
import PyPlot
mpl = plt.matplotlib
# mpl.style.use(GRASS.moddir * "fig.mplstyle")
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

# get command line args and output directories
# include(joinpath(abspath(@__DIR__), "paths.jl"))
plotdir = string(abspath(joinpath("figures", "NEID_sim_comparison")))
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
#line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]

#determine timestamps
# # October Eclipse
# path = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
# timestamps_full = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
# timestamps = timestamps_full[16:length(timestamps_full)-150]
#April Eclipse
path = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2024/04/08/"
timestamps_full = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_April/data/neid_april_data.csv", DataFrame)[!, "filename"]
timestamps = timestamps_full[101:length(timestamps_full)-10]

# read in optimized depths
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

df_tuned = CSV.read(joinpath(datadir, "tuned_params.csv"), DataFrame)

full_pixels = 2048:7168
# wavelength of line to synthesize/compare to neid
for i in 8:length(files)
    println(">>> Running " * line_names[i] * "...")

    # get properties from line
    line_name = line_names[i]
    airwav = lp.λrest[i]
    airvac = airwav .* refractive_index
    depth = lp.depth[i]

    order_index = orders[i]
    min_wav = airvac - 1.5
    max_wav = airvac + 1.5

   # get depth from neid spectrum
    buff = 0.2
    if contains("FeI_5383", line_name)
        buff = 0.3
    elseif contains("FeI_5434", line_name)
        buff = 0.3
    elseif contains("FeI_5436.3", line_name)
        buff = 0.135
    elseif contains("FeI_5382", line_name)
        buff = 0.2
    elseif contains("FeI_5576", line_name)
        buff = 0.25
    elseif contains("CaI_6169.0", line_name)
        buff = 0.25
    elseif contains("FeI_6301", line_name)
        buff = 0.25
    end

    # get the depth for the simulation
    sim_depth = df[i, "optimized_depth"]

    # simulate the spectrum
    lines = [airwav]
    depths = [sim_depth]
    templates = [files[i]]
    resolution = 11e4
    spec = SpecParams(lines=lines, depths=depths, templates=templates,
                      resolution=resolution, buffer = 1.5, oversampling=4.0)
    N = 5
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=5)

    # simulate the spectrum
    wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, obs_long, obs_lat, alt, lines ./ 10.0, time_stamps, "Optical", verbose=true, use_gpu=false)

    #convert from airwav to vacuum wavelength 
    spec.lambdas .= spec.lambdas * refractive_index
    spec.lines .= spec.lines * refractive_index

    for j in 1:1#length(timestamps)
        spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)

        #find where line is
        chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
        #pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
        pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
        wavs_neid = chunk.λ[pixels]
        flux_neid = chunk.flux[pixels]
        flux_neid ./= maximum(flux_neid)

        idxl = findfirst(x -> x .>= airvac - buff, wavs_neid)
        idxr = findfirst(x -> x .>= airvac + buff, wavs_neid)
        neid_bot = minimum(view(flux_neid, idxl:idxr))
        neid_depth = 1.0 - neid_bot

        # interpolate neid on synth wavelength grid
        itp = GRASS.linear_interp(wavs_neid, flux_neid, bc=NaN)
        flux_neid = itp.(wavs_sim)
        wavs_neid = copy(wavs_sim)

        flux_sim = flux_sim[((j-1)*length(wavs_sim) + 1):(length(wavs_sim)*j)]

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
    
        v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim, spec)

        v_grid_neid, ccf_neid = GRASS.calc_ccf(wavs_neid, flux_neid, spec)

        plt.figure()
        plt.plot(v_grid_neid, ccf_neid, label = "NEID")
        plt.plot(v_grid_sim, ccf_sim, label = "GRASS")
        plt.legend()
        plt.savefig("ccf_test_$i.png")

        # deal with annoying line blend
        if contains("FeI_6302", line_name)
            amin_ccf = argmin(ccf_neid)
            idx_b = findfirst(x -> x .> 0.925, ccf_neid[amin_ccf:end]) + amin_ccf
            ccf_neid = ccf_neid[1:idx_b]
            v_grid_neid = v_grid_neid[1:idx_b]
        end

        # get bisectors
        top = 0.9
        vel_sim, int_sim = GRASS.calc_bisector(v_grid_sim, ccf_sim, nflux=50, top=top)
        vel_neid, int_neid = GRASS.calc_bisector(v_grid_neid, ccf_neid, nflux=50, top=top)

        # compute velocity as mean bisector between N and M % depth
        N = 0.20
        M = 0.70
        idx1 = findfirst(x -> x .>= N * sim_depth + minimum(flux_sim), int_sim)
        idx2 = findfirst(x -> x .>= M * sim_depth + minimum(flux_sim), int_sim)
        if isnothing(idx2)
            idx2 = findfirst(x -> x .>= 0.9, int_sim)
        end
        rv_sim = mean(view(vel_sim, idx1:idx2))

        idx1 = findfirst(x -> x .>= N * neid_depth + neid_bot, int_neid)
        idx2 = findfirst(x -> x .>= M * neid_depth + neid_bot, int_neid)
        if isnothing(idx2)
            idx2 = findfirst(x -> x .>= 0.9, int_neid)
        end
        rv_neid = mean(view(vel_neid, idx1:idx2))

        # transform to lab frame
        vel_neid .-= rv_neid
        vel_sim .-= rv_sim
        wavs_neid ./= calc_doppler_factor(rv_neid)
        wavs_sim ./= calc_doppler_factor(rv_sim)

        # interpolate neid onto synthetic wavelength grid
        itp = GRASS.linear_interp(wavs_neid, flux_neid)
        flux_neid = itp.(wavs_sim)
        wavs_neid = copy(wavs_sim)

        # recompute bisectors b/c of interpolation
        v_grid_neid, ccf_neid = GRASS.calc_ccf(wavs_neid, flux_neid, spec)

        v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim, spec)

        if contains("FeI_6302", line_name)
            amin_ccf = argmin(ccf_neid)
            idx_b = findfirst(x -> x .> 0.925, ccf_neid[amin_ccf:end]) + amin_ccf
            ccf_neid = ccf_neid[1:idx_b]
            v_grid_neid = v_grid_neid[1:idx_b]
        end

        # get bisectors
        vel_neid, int_neid = GRASS.calc_bisector(v_grid_neid, ccf_neid, nflux=50, top=top)
        vel_sim, int_sim = GRASS.calc_bisector(v_grid_sim, ccf_sim, nflux=50, top=top)

        # find mean velocities in order to align bisectors
        idx1 = findfirst(x -> x .>= N * sim_depth + minimum(flux_sim), int_sim)
        idx2 = findfirst(x -> x .>= M * sim_depth + minimum(flux_sim), int_sim)
        if isnothing(idx2)
            idx2 = findfirst(x -> x .>= 0.9, int_sim)
        end
        rv_sim = mean(view(vel_sim, idx1:idx2))

        idx1 = findfirst(x -> x .>= N * neid_depth + neid_bot, int_neid)
        idx2 = findfirst(x -> x .>= M * neid_depth + neid_bot, int_neid)
        if isnothing(idx2)
            idx2 = findfirst(x -> x .>= 0.9, int_neid)
        end
        rv_neid = mean(view(vel_neid, idx1:idx2))

        # align the bisectors
        vel_sim .-= rv_sim
        vel_neid .-= rv_neid

        # get the bisector residuals
        bis_resids = vel_neid .- vel_sim
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

        plt.figure()
        plt.plot(wavs_neid, flux_neid, color = "b")
        plt.scatter(wavs_neid, flux_neid, label = "NEID", color = "b", s = 5)
        plt.plot(wavs_sim, flux_sim, color = "r")
        plt.scatter(wavs_sim, flux_sim, label = "GRASS", color = "r", s = 5)
        plt.legend()
        # plt.axvline(x=airvac)
        plt.axvline(x=airvac - buff)
        plt.axvline(x=airvac + buff)
        plt.title(line_names[i])
        plt.gca().invert_xaxis()
        plt.savefig("spectra_zoom_test_$i.png")

        # # big function for plotting
        # function comparison_plots()
        #     # make plot objects
        #     fig = plt.figure(figsize=(6.4,4.8))
        #     gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
        #     ax1 = fig.add_subplot(gs[1])
        #     ax2 = fig.add_subplot(gs[2])

        #     # plot spectra
        #     ax1.plot(wavs_sim, flux_sim, marker="o", c="black", ms=3.0, lw=1.5, markevery=2, label="Synthetic")
        #     ax1.plot(wavs_neid, flux_neid, marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=2, label="neid")
        #     # ax1.plot(wavs_neid, flux_mod, alpha=0.9, marker="o", c=colors[2], ms=2.0, lw=1.0, markevery=2, label=L"{\rm Cleaned\ neid}")

        #     # plot resids
        #     ax2.plot(wavs_sim, flux_neid .- flux_sim, c=colors[1], marker="s", ms=2.0, lw=0, markevery=2)

        #     # find limits
        #     idx_min = argmin(flux_sim)
        #     idx1 = idx_min - findfirst(x -> x .> 0.95, flux_sim[idx_min:-1:1])
        #     idx2 = idx_min + findfirst(x -> x .> 0.95, flux_sim[idx_min:end])

        #     # set limits
        #     min_idx = argmin(flux_neid)
        #     ax1.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
        #     ax1.set_ylim(minimum(flux_sim) - 0.1, 1.1)
        #     ax2.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
        #     ax2.set_ylim(-0.0575, 0.0575)

        #     # set tick labels, axis labels, etc.
        #     ax1.set_xticklabels([])
        #     ax1.set_ylabel(L"{\rm Normalized\ Flux}", labelpad=15)
        #     ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
        #     ax2.set_ylabel(L"{\rm neid\ -\ GRASS}")
        #     ax1.legend()
        #     fig.tight_layout()

        #     # set the title
        #     title = replace(line_name, "_" => "\\ ")
        #     idx = findfirst('I', title)
        #     title = title[1:idx-1] * "\\ " * title[idx:end] * "\\ \\AA"
        #     ax1.set_title(("\${\\rm " * title * "}\$"))

        #     # save the plot
        #     fig.subplots_adjust(wspace=0.05)
        #     fig.savefig(joinpath(plotdir, line_name * "_line.png"))
        #     plt.clf(); plt.close()

        #     # plot the bisectors
        #     fig = plt.figure(figsize=(6.4,4.8))
        #     gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
        #     ax1 = fig.add_subplot(gs[1])
        #     ax2 = fig.add_subplot(gs[2])

        #     # plot bisectors
        #     ax1.plot(vel_sim[4:end], int_sim[4:end], marker="o", color="black", ms=3.0, lw=2.0, markevery=1, label="Synthetic")
        #     ax1.plot(vel_neid[4:end], int_neid[4:end], marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=1, label="neid")

        #     # plot residuals
        #     ax2.plot(vel_neid[4:end] .- vel_sim[4:end], int_neid[4:end], c=colors[1], marker="s", ms=2.0, lw=0.0, markevery=1)
        #     # set tick labels, axis labels, etc.
        #     ax2.set_yticklabels([])
        #     ax2.yaxis.tick_right()
        #     ax1.set_ylim(minimum(flux_sim) - 0.05, 1.05)
        #     ax2.set_xlim(-35, 35)
        #     ax2.set_ylim(minimum(flux_sim) - 0.05, 1.05)
        #     ax1.set_xlabel(L"{\rm Relative\ Velocity\ (m\ s^{-1})}", fontsize=13)
        #     ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
        #     ax2.set_xlabel(L"{\rm neid\ -\ GRASS\ (m\ s^{-1})}", fontsize=13)
        #     ax1.legend(labelspacing=0.25)
        #     ax2.legend(loc="upper right", prop=Dict("size"=>12.5), labelspacing=0.25)

        #     # set the title
        #     fig.suptitle(("\${\\rm " * title * "}\$"), y=0.95)

        #     # save the plot
        #     fig.subplots_adjust(hspace=0.05)
        #     fig.savefig(joinpath(plotdir, line_name * "_bisector.png"))
        #     plt.clf(); plt.close()
        #     return nothing
        # end

        # comparison_plots()
    end
end
