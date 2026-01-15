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
using EchelleCCFs: λ_air_to_vac, λ_vac_to_air, calc_doppler_factor, MeasureRvFromCCFQuadratic as QuadraticFit
using RvSpectMLBase
using EchelleInstruments
# plotting
using LaTeXStrings
import PyPlot
plt = PyPlot
mpl = plt.matplotlib
import PyPlot: matplotlib
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

GRASS.get_kernels()

using PyCall
py"""
import sys
sys.path.append('.')
"""
neid_lsf = pyimport("NEID_LSF")
neid_asymmetric_lsf = pyimport("NeidLsf")
np = pyimport("numpy")
scipy_interp = pyimport("scipy.interpolate")
CubicSpline = scipy_interp.CubicSpline
mdates = pyimport("matplotlib.dates")

#determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)
λrest = GRASS.get_rest_wavelength(lp)
depth = GRASS.get_depth(lp)
lfile = GRASS.get_file(lp)

#line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
variable_names = ["lambda_min", "lambda_max", "flux_min", "pixel_mean", "pixels", "neid_wavelength"]
ext_coeff_array = [0.15452995224327976, 0.15256098077094832, 0.14720055859068512, 0.154895798933504, 0.15181381895180662, 0.15107508233588227, 0.15116772762156633, 0.14882114581650618, 0.14865707189399568, 0.1494903120065096, 0.16011027092744037, 0.15593033972594958, 0.14195968590211427, 0.15401904166429853, 0.1277699772941639, 0.12709315507233226, 0.12820346527304866, 0.11702310600015708, 0.1435320747844216, 0.12380490304619193, 0.12450734135297492, 0.12101777355247835]

datadir = string(abspath("data"))
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938

function neid_bis(neid_timestamps, timestamps, path, LD_type)
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    N = 197
    Nt = length(neid_timestamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=40)

    full_pixels = 2048:7168
    # wavelength of line to synthesize/compare to neid
    for i in 1:length(line_names)
        # get properties from line
        line_name = line_names[i]
        order_index = orders[i]
        min_wav = vacwav[i] - 5.0
        max_wav = vacwav[i] + 5.0

        buff = 0.3
        if line_name == "FeI_5383"
            buff = 0.55
        end
        if line_name == "FeI_5436.3"
            buff = 0.25
        end

        #simulate the spectrum
        lines = [λrest[i]]
        templates = [lfile[i]]
        depths = [df[i, "optimized_depth"]]
        variability = trues(length(lines))  # whether or not the bisectors should "dance"
        blueshifts = zeros(length(lines))   # set convective blueshift value
        resolution = 7e5

        data = jldopen("data/NEID_convolution_info$(i).jld2", "r") do file
            Dict(var => read(file, var) for var in variable_names)
        end
        pixel_mean = data["pixel_mean"]

        # make the spec composite type instances
        spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
                        blueshifts=blueshifts, templates=templates, resolution=resolution) 

        # simulate the spectrum 
        wavs_sim, flux_sim = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, obs_long, obs_lat, alt, time_stamps, lines, ext_coeff_array[i], true, true)
        #convert from airwav to vacuum wavelength 
        wavs_sim .= λ_air_to_vac.(wavs_sim)

        for j in 1:length(time_stamps)
            outspec_t = flux_sim[:, j]
            outspec_t = outspec_t ./ maximum(outspec_t)                                                            

            model_pixels, flux_sim_j, dlam_dpix = neid_asymmetric_lsf.lsf_model_export(orders[i], round(Int, pixel_mean[j]), wavs_sim, outspec_t, vacwav[i])
            wavs_sim = (dlam_dpix .* model_pixels) .+ vacwav[i]
            wavs_sim .= λ_vac_to_air.(wavs_sim)

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
            pixel_normal = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            pixels_inner = NEID.find_pixels_for_line_in_chunk(chunk, vacwav[i] - buff, vacwav[i] + buff)
            wavs_neid = chunk.λ[pixels_inner]
            wavs_neid .= λ_vac_to_air.(wavs_neid)
            chunk_flux = chunk.flux[pixels_inner]
            flux_neid = chunk_flux ./ maximum(chunk.flux[pixel_normal])

            neid_bot = minimum(flux_neid)
            neid_depth = 1.0 - neid_bot

            sim_bot = minimum(flux_sim_j)
            sim_depth = 1.0 - sim_bot
    
            # interpolate neid on synth wavelength grid
            itp = GRASS.linear_interp(wavs_neid, flux_neid, bc=NaN)
            flux_neid = itp.(wavs_sim)
            wavs_neid = copy(wavs_sim)

            v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim_j, lines, [sim_depth], 11e4)
            v_grid_neid, ccf_neid = GRASS.calc_ccf(wavs_neid, flux_neid, lines, [neid_depth], 11e4)

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
            idx1 = findfirst(x -> x .>= N * sim_depth + minimum(flux_sim_j), int_sim)
            idx2 = findfirst(x -> x .>= M * sim_depth + minimum(flux_sim_j), int_sim)
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
            v_grid_neid, ccf_neid = GRASS.calc_ccf(wavs_neid, flux_neid, lines, [neid_depth], 11e4)
            v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim_j, lines, [sim_depth], 11e4)

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
            idx1 = findfirst(x -> x .>= N * sim_depth + minimum(flux_sim_j), int_sim)
            idx2 = findfirst(x -> x .>= M * sim_depth + minimum(flux_sim_j), int_sim)
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

            # big function for plotting
            function comparison_plots()
                # make plot objects
                fig = plt.figure(figsize=(6.4,4.8))
                gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
                ax1 = fig.add_subplot(gs[1])
                ax2 = fig.add_subplot(gs[2])
            
                # plot spectra
                ax1.plot(wavs_sim, flux_sim_j, marker="o", c="black", ms=3.0, lw=2, markevery=15, label="Synthetic")
                ax1.plot(wavs_neid, flux_neid, marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=15, label="NEID")
                # plot resids
                ax2.plot(wavs_sim, flux_neid .- flux_sim_j, c=colors[1], marker="s", ms=2.0, lw=0, markevery=15)
            
                # find limits
                idx_min = argmin(flux_sim_j)
                idx1 = idx_min - findfirst(x -> x .> 0.95, flux_sim_j[idx_min:-1:1])
                idx2 = idx_min + findfirst(x -> x .> 0.95, flux_sim_j[idx_min:end])
            
                # set limits
                min_idx = argmin(flux_neid)
                ax1.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
                ax1.set_ylim(minimum(flux_sim_j) - 0.1, 1.1)
                ax2.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
                ax2.set_ylim(-0.0275, 0.0275)
            
                # set tick labels, axis labels, etc.
                ax1.set_xticklabels([])
                fmt = matplotlib.ticker.ScalarFormatter()
                fmt.set_useOffset(false)
                fmt.set_scientific(false)

                ax2.xaxis.set_major_formatter(fmt)   # lower plot
                mticker = PyCall.pyimport("matplotlib.ticker")
                ax2.xaxis.set_major_locator(mticker.MaxNLocator(5))

                ax1.set_ylabel(L"{\rm Normalized\ Intensity}", labelpad=15, fontsize=13)
                ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}", fontsize=13)
                ax2.set_ylabel(L"{\rm NEID\ -\ GRASS}", fontsize=13)
                ax1.tick_params(axis="both", labelsize=13)
                ax2.tick_params(axis="both", labelsize=13)
                ax1.legend(fontsize=13)
                fig.tight_layout()
            
                # set the title
                title = replace(line_name, "_" => "\\ ")
                idx = findfirst('I', title)
                title = title[1:idx-1] * "\\ " * title[idx:end] * "\\ \\AA"
                ax1.set_title(("\${\\rm " * title * "}\$"), fontsize=13)
            
                # save the plot
                fig.subplots_adjust(wspace=0.05)
                fig.savefig("line_profiles/f10_$(i).pdf", bbox_inches="tight")
                plt.clf(); plt.close()
            
                # plot the bisectors
                fig = plt.figure(figsize=(6.4,4.8))
                gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
                ax1 = fig.add_subplot(gs[1])
                ax2 = fig.add_subplot(gs[2])
            
                # plot bisectors
                ax1.plot(vel_sim[4:end], int_sim[4:end], marker="o", color="black", ms=3.0, lw=2.0, markevery=1, label="Synthetic")
                ax1.plot(vel_neid[4:end], int_neid[4:end], marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=1, label="NEID")
            
                # plot residuals
                ax2.plot(vel_neid[4:end] .- vel_sim[4:end], int_neid[4:end], c=colors[1], marker="s", ms=2.0, lw=0.0, markevery=1)
                # set tick labels, axis labels, etc.
                ax2.set_yticklabels([])
                ax2.yaxis.tick_right()
                ax1.set_ylim(minimum(flux_sim_j) - 0.05, 1.05)
                ax2.set_xlim(-35, 35)
                ax2.set_ylim(minimum(flux_sim_j) - 0.05, 1.05)
                ax1.tick_params(axis="both", labelsize=13)
                ax2.tick_params(axis="both", labelsize=13)
                ax1.set_xlabel(L"{\rm Relative\ Velocity\ (m/s)}", fontsize=13)
                ax1.set_ylabel(L"{\rm Normalized\ Intensity}", fontsize=13)
                ax2.set_xlabel(L"{\rm NEID\ -\ GRASS\ (m/s)}", fontsize=13)
                ax1.legend(fontsize=13)
            
                # set the title
                fig.suptitle(("\${\\rm " * title * "}\$"), y=0.95, fontsize=13)
            
                # save the plot
                fig.subplots_adjust(hspace=0.05)
                fig.savefig("line_profiles/f10_$(i + 22).pdf", bbox_inches="tight")
                plt.clf(); plt.close()
                return nothing
            end

            comparison_plots()
        end
    end
end

neid_timestamps = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
timestamps = timestamps_full_october[16:length(timestamps_full_october)-150]

neid_bis(neid_timestamps[130:130], timestamps[130:130], path_october, "SSD_4parameter")
