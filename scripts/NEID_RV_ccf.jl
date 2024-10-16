using CSV
using JLD2
using GRASS
using SPICE
using FITSIO
using FileIO
using Statistics
using DataFrames
using EchelleCCFs
using RvSpectMLBase
using EchelleInstruments
using EchelleCCFs: λ_air_to_vac, λ_vac_to_air
# plotting
using LaTeXStrings
import PyPlot
plt = PyPlot
mpl = plt.matplotlib
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

#determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)

#line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]

datadir = string(abspath("data"))
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

neid_timestamps_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
# October Eclipse
path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
timestamps_october = timestamps_full_october[16:length(timestamps_full_october)-150]

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

function read_iag_atlas(;isolate::Bool=true, airwav::Float64=5434.5232, buffer::Float64=1.0)
    file = datadir * "/spvis.dat"

    # read in the IAG atlas
    iag = CSV.read(file, DataFrame, ignorerepeated=true, delim="|", skipto=5,
                   footerskip=1, header=["wavenum", "nflux", "flux"])
    # convert wavenumber to wavelength in angstroms
    wavs = (1.0 ./ iag.wavenum) * 1e8
    # reverse to deal with conversion of units
    reverse!(wavs)
    reverse!(iag.nflux)

    # isolate region around line
    if isolate
        vacwav = λ_air_to_vac(airwav)
        ind1 = findfirst(x -> x .> vacwav-buffer, wavs)
        ind2 = findfirst(x -> x .> vacwav+buffer, wavs[ind1:end]) + ind1
        return (view(wavs,ind1:ind2)), view(iag.nflux, ind1:ind2)
    else
        return wavs, iag.nflux
    end
end

function last_timestamp_lines(line_names, airwav, vacwav, orders, neid_timestamps, timestamps, path, filename, LD_type)
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    variable_names = ["zenith_mean", "dA_total_proj", "idx1", "idx3", "mu_grid", "z_rot_sub", "mu", "ax_codes", "dA", "N"]

    # Open the JLD2 file and read the variables into a dictionary
    data = jldopen("data/solar_disk/$(filename).jld2", "r") do file
        Dict(var => read(file, var) for var in variable_names)
    end

    zenith_mean = deepcopy(data["zenith_mean"])
    dA_total_proj = deepcopy(data["dA_total_proj"])
    idx1 = deepcopy(data["idx1"])
    idx3 = deepcopy(data["idx3"])
    mu_grid = deepcopy(data["mu_grid"])
    z_rot_sub = deepcopy(data["z_rot_sub"])
    stored_μs = deepcopy(data["mu"])
    stored_ax_codes = deepcopy(data["ax_codes"])
    stored_dA = deepcopy(data["dA"])

    N = data["N"]
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    full_pixels = 2048:7168
    for i in 1:length(line_names)
        order_index = orders[i]
        min_wav = vacwav[i] - 2
        max_wav = vacwav[i] + 2

        #get the depth for the simulation
        sim_depth = df[i, "optimized_depth"]
        #simulate the spectrum
        lines = [airwav[i]]
        depths = [sim_depth]
        templates = [line_names[i]]
        resolution = 7e5
        spec = GRASS.SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution, oversampling=4.0)
    
        # simulate the spectrum 
        wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_type, zenith_mean,
                                                            dA_total_proj, idx1, idx3, mu_grid, z_rot_sub,
                                                            stored_μs, stored_ax_codes, stored_dA, [0.0], ext_toggle = false, verbose=true, use_gpu=false)
        
        # convolve GRASS spectrum to NEID resolution
        wavs_sim, flux_sim = GRASS.convolve_gauss(wavs_sim, flux_sim, new_res=11e4, oversampling=4.0)
        wavs_sim .= λ_air_to_vac.(wavs_sim)

        for j in length(timestamps):length(timestamps)
            #GRASS flux for given timestamps
            flux_sim_j = flux_sim[((j-1)*length(wavs_sim) + 1):(length(wavs_sim)*j)]

            #IAG flux 
            wavs_iag0, flux_iag0 = read_iag_atlas(isolate=true, airwav=airwav[i], buffer=2.0)
            flux_iag0 .-= 1 - maximum(flux_sim_j)

            f = FITS(joinpath(path, timestamps[j]))
            header = read_header(f[1])
            bc_ms = (header["SSBRV0$order_index"]) * 1000

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            #find where line is - NEID flux
            chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
            pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            chunk_flux_full_o = chunk.flux[pixels]
            chunk_flux_full = chunk_flux_full_o ./ maximum(chunk_flux_full_o)
            chunk_flux_full .-= 1 - maximum(flux_sim_j)
            chunck_vac_wav = chunk.λ[pixels] ./ ((-bc_ms)/GRASS.c_ms + 1)

            #plot spectrum
            plt.figure()
            plt.plot(wavs_iag0, flux_iag0, label = "IAG", color = "g")
            plt.plot(chunck_vac_wav, chunk_flux_full, label = "NEID", color = "b")
            plt.scatter(chunck_vac_wav, chunk_flux_full, color = "b", s = 5)
            plt.plot(wavs_sim, flux_sim_j, label = "GRASS", color = "r")
            plt.scatter(wavs_sim, flux_sim_j, color = "r", s = 5)
            plt.axvline(x = vacwav[i])
            plt.legend()

            pixels = NEID.find_pixels_for_line_in_chunk(chunk, vacwav[i] - 0.25, vacwav[i] + 0.25)
            chunk_flux = chunk.flux[pixels]
            chunk_flux ./= maximum(chunk_flux_full_o)
            chunk_flux .-= 1 - maximum(flux_sim_j)
            chunck_vac_wav = chunk.λ[pixels] ./ ((-bc_ms)/GRASS.c_ms + 1)

            sim_shift_wl = (round(wavs_sim[findmin(flux_sim_j)[2]] - vacwav[i]; digits = 6))
            data_shift_wl = (round(chunck_vac_wav[findmin(chunk_flux)[2]] - vacwav[i]; digits = 6))
            sim_shift_vel = round(sim_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)
            data_shift_vel = round(data_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)

            plt.text(wavs_sim[120], 0.55, "Sim. Shift $sim_shift_wl Å")
            plt.text(wavs_sim[120], 0.45, "Sim. Shift $sim_shift_vel m/s")
            plt.text(wavs_sim[120], 0.5, "Data Shift $data_shift_wl Å")
            plt.text(wavs_sim[120], 0.4, "Data Shift $data_shift_vel m/s")

            plt.title("air wav: $(line_names[i])")
            plt.gca().invert_xaxis()
            plt.savefig("eclipse_figures/Spectrum/October_zoom_out_comp/spectra_$(i).png")
        end
    end
end

function GRASS_comparison(line_names, airwav, vacwav, orders, neid_timestamps, timestamps, path, filename, LD_type)
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    variable_names = ["zenith", "dA_total_proj", "idx1", "idx3", "mu_grid", "z_rot_sub", "mu", "ax_codes", "dA", "N"]

    # Open the JLD2 file and read the variables into a dictionary
    data = jldopen("data/solar_disk/$(filename).jld2", "r") do file
        Dict(var => read(file, var) for var in variable_names)
    end

    N = data["N"]
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    full_pixels = 2048:7168
    #iterate through lines and determine line RV for eclipse EM curve
    for i in 1:length(line_names)
        if i == 11
            continue
        end

        zenith_mean = deepcopy(data["zenith"])
        dA_total_proj = deepcopy(data["dA_total_proj"])
        idx1 = deepcopy(data["idx1"])
        idx3 = deepcopy(data["idx3"])
        mu_grid = deepcopy(data["mu_grid"])
        z_rot_sub = deepcopy(data["z_rot_sub"])
        stored_μs = deepcopy(data["mu"])
        stored_ax_codes = deepcopy(data["ax_codes"])
        stored_dA = deepcopy(data["dA"])

        order_index = orders[i]
        min_wav = vacwav[i] - 2
        max_wav = vacwav[i] + 2
        line_name = line_names[i]

        #get the depth for the simulation
        sim_depth = df[i, "optimized_depth"]
        #simulate the spectrum
        lines = [airwav[i]]
        depths = [sim_depth]
        templates = [line_names[i]]
        resolution = 7e5
        spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution, oversampling=4.0)
    
        # simulate the spectrum 
        wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_type, zenith_mean,
                                                            dA_total_proj, deepcopy(idx1), idx3, mu_grid, z_rot_sub,
                                                            stored_μs, stored_ax_codes, stored_dA, [0.0], ext_toggle = false, verbose=true, use_gpu=false)
        # convolve GRASS spectrum to NEID resolution
        wavs_sim, flux_sim = GRASS.convolve_gauss(wavs_sim, flux_sim, new_res=11e4, oversampling=4.0)
        wavs_sim .= λ_air_to_vac.(wavs_sim)

        #convert from airwav to vacuum wavelength 
        spec.lambdas .= λ_air_to_vac.(spec.lambdas)
        spec.lines .= λ_air_to_vac.(spec.lines)

        for j in length(timestamps):length(timestamps)
            #GRASS flux for given timestamps
            flux_sim_j = flux_sim[((j-1)*length(wavs_sim) + 1):(length(wavs_sim)*j)]

            f = FITS(joinpath(path, timestamps[j]))
            header = read_header(f[1])
            bc_ms = (header["SSBRV0$order_index"]) * 1000

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            #find where line is - NEID flux
            chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
            pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            chunk_flux_full_o = chunk.flux[pixels]
            chunk_flux_full = chunk_flux_full_o ./ maximum(chunk_flux_full_o)
            chunk_flux_full .-= 1 - maximum(flux_sim_j)
            chunck_vac_wav = chunk.λ[pixels] ./ ((-bc_ms)/GRASS.c_ms + 1)

            # interpolate neid on synth wavelength grid
            itp = GRASS.linear_interp(chunck_vac_wav, chunk_flux_full, bc=NaN)
            chunk_flux_full = itp.(wavs_sim)
            chunck_vac_wav = copy(wavs_sim)

            v_grid_cpu_neid, ccf_cpu_neid = GRASS.calc_ccf(chunck_vac_wav, chunk_flux_full, lines, [maximum(chunk_flux_full) - minimum(chunk_flux_full)], 11e4)
            v_grid_cpu_sim, ccf_cpu_sim = GRASS.calc_ccf(wavs_sim, flux_sim_j, spec)

            # deal with annoying line blend
            if contains("FeI_6302", line_name)
                amin_ccf = argmin(ccf_cpu_neid)
                idx_b = findfirst(x -> x .> 0.925, ccf_cpu_neid[amin_ccf:end]) + amin_ccf
                ccf_cpu_neid = ccf_cpu_neid[1:idx_b]
                v_grid_cpu_neid = v_grid_cpu_neid[1:idx_b]
            end

            # get bisectors
            top = 0.9
            vel_sim, int_sim = GRASS.calc_bisector(v_grid_cpu_sim, ccf_cpu_sim, nflux=50, top=top)
            vel_neid, int_neid = GRASS.calc_bisector(v_grid_cpu_neid, ccf_cpu_neid, nflux=50, top=top)

            # big function for plotting
            function comparison_plots()
                # make plot objects
                fig = plt.figure(figsize=(6.4,4.8))
                gs = mpl.gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1], figure=fig, hspace=0.05)
                ax1 = fig.add_subplot(gs[1])
                ax2 = fig.add_subplot(gs[2])
            
                # plot spectra
                ax1.plot(wavs_sim, flux_sim_j, marker="o", c="black", ms=3.0, lw=1.5, markevery=2, label="Synthetic")
                ax1.plot(chunck_vac_wav, chunk_flux_full, marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=2, label="neid")
                # plot resids
                ax2.plot(wavs_sim, chunk_flux_full .- flux_sim_j, c=colors[1], marker="s", ms=2.0, lw=0, markevery=2)
            
                # find limits
                idx_min = argmin(flux_sim_j)
                idx1 = idx_min - findfirst(x -> x .> 0.95, flux_sim_j[idx_min:-1:1])
                idx2 = idx_min + findfirst(x -> x .> 0.95, flux_sim_j[idx_min:end])
            
                # set limits
                min_idx = argmin(chunk_flux_full)
                ax1.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
                ax1.set_ylim(minimum(flux_sim_j) - 0.1, 1.1)
                ax2.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
                ax2.set_ylim(-0.0575, 0.0575)
            
                # set tick labels, axis labels, etc.
                ax1.set_xticklabels([])
                ax1.set_ylabel(L"{\rm Normalized\ Flux}", labelpad=15)
                ax2.set_xlabel(L"{\rm Wavelength\ (\AA)}")
                ax2.set_ylabel(L"{\rm neid\ -\ GRASS}")
                ax1.legend()
                fig.tight_layout()
            
                # set the title
                title = replace(line_name, "_" => "\\ ")
                idx = findfirst('I', title)
                title = title[1:idx-1] * "\\ " * title[idx:end] * "\\ \\AA"
                ax1.set_title(("\${\\rm " * title * "}\$"))

                ax1.axvline(x = vacwav[i])
                sim_shift_wl = (round(wavs_sim[findmin(flux_sim_j)[2]] - vacwav[i]; digits = 6))
                data_shift_wl = (round(chunck_vac_wav[idx_min-10:idx_min+10][findmin(chunk_flux_full[idx_min-10:idx_min+10])[2]] - vacwav[i]; digits = 6))
                sim_shift_vel = round(sim_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)
                data_shift_vel = round(data_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)
    
                ax1.text(wavs_sim[120], 0.75, "Sim. Shift $sim_shift_wl Å")
                ax1.text(wavs_sim[120], 0.70, "Sim. Shift $sim_shift_vel m/s")
                ax1.text(wavs_sim[120], 0.65, "Data Shift $data_shift_wl Å")
                ax1.text(wavs_sim[120], 0.60, "Data Shift $data_shift_vel m/s")
            
                # save the plot
                fig.subplots_adjust(wspace=0.05)
                fig.savefig(joinpath("eclipse_figures/Spectrum/October_zoom_in_comp", line_name * "_line.png"))
                plt.clf(); plt.close()
            
                # plot the bisectors
                fig = plt.figure(figsize=(6.4,4.8))
                gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
                ax1 = fig.add_subplot(gs[1])
                ax2 = fig.add_subplot(gs[2])
            
                # plot bisectors
                ax1.plot(vel_sim[4:end], int_sim[4:end], marker="o", color="black", ms=3.0, lw=2.0, markevery=1, label="Synthetic")
                ax1.plot(vel_neid[4:end], int_neid[4:end], marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=1, label="neid")
            
                # plot residuals
                ax2.plot(vel_neid[4:end] .- vel_sim[4:end], int_neid[4:end], c=colors[1], marker="s", ms=2.0, lw=0.0, markevery=1)
                # set tick labels, axis labels, etc.
                ax2.set_yticklabels([])
                ax2.yaxis.tick_right()
                ax1.set_ylim(minimum(flux_sim_j) - 0.05, 1.05)
                ax2.set_xlim(-35, 35)
                ax2.set_ylim(minimum(flux_sim_j) - 0.05, 1.05)
                ax1.set_xlabel(L"{\rm Relative\ Velocity\ (m\ s^{-1})}", fontsize=13)
                ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
                ax2.set_xlabel(L"{\rm neid\ -\ GRASS\ (m\ s^{-1})}", fontsize=13)
                ax1.legend(labelspacing=0.25)
                ax2.legend(loc="upper right", prop=Dict("size"=>12.5), labelspacing=0.25)
            
                # set the title
                fig.suptitle(("\${\\rm " * title * "}\$"), y=0.95)
            
                # save the plot
                fig.subplots_adjust(hspace=0.05)
                fig.savefig(joinpath("eclipse_figures/Spectrum/October_zoom_in_comp", line_name * "_bisector.png"))
                plt.clf(); plt.close()
                return nothing
            end
            comparison_plots()
        end
    end
end

function line_rvs_ccf(line_names, vacwav, orders, timestamps, path)

    full_pixels = 2048:7168

    resolution = 11e4

    RV_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)
    RV_error_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)
    #iterate through lines and determine line RV for eclipse EM curve
    Threads.@threads for i in 1:length(line_names)
        order_index = orders[i]
        min_wav = vacwav[i] - 2
        max_wav = vacwav[i] + 2

        sim_depth = df[i, "optimized_depth"]
        lines = [vacwav[i]]

        RV_list = Vector{Float64}(undef,length(timestamps)...)
        RV_error_list = Vector{Float64}(undef,length(timestamps)...)
        Threads.@threads for j in 1:length(timestamps)
            f = FITS(joinpath(path, timestamps[j]))
            header = read_header(f[1])
            bc_ms = (header["SSBRV0$order_index"]) * 1000

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            #find where line is - NEID flux
            chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
            pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            chunk_flux_full_o = chunk.flux[pixels]
            chunk_flux_full = chunk_flux_full_o ./ maximum(chunk_flux_full_o)
            chunck_vac_wav = chunk.λ[pixels]

            v_grid_cpu, ccf_cpu = GRASS.calc_ccf(chunck_vac_wav, chunk_flux_full, lines, [maximum(chunk_flux_full) - minimum(chunk_flux_full)], resolution)
            rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

            RV_list[j] = rvs_cpu
            RV_error_list[j] = sigs_cpu

            if line_names[i] == "FeI_5383"
                fig = plt.figure(figsize=(6.4,4.8))
                gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
                ax1 = fig.add_subplot(gs[1])
                ax2 = fig.add_subplot(gs[2])

                ax1.plot(chunck_vac_wav ./ ((-bc_ms)/GRASS.c_ms + 1), chunk_flux_full, c=colors[1])
                ax1.scatter(chunck_vac_wav ./ ((-bc_ms)/GRASS.c_ms + 1), chunk_flux_full, c=colors[1])
                ax1.axvline(x = vacwav[i], c=colors[1])
                # find limits
                idx_min = argmin(chunk_flux_full)
                idx1 = idx_min - findfirst(x -> x .> 0.95, chunk_flux_full[idx_min:-1:1])
                idx2 = idx_min + findfirst(x -> x .> 0.95, chunk_flux_full[idx_min:end])
                # set limits
                ax1.set_xlim(chunck_vac_wav[idx1-20], chunck_vac_wav[idx2+20])
                ax1.set_ylim(minimum(chunk_flux_full) - 0.1, 1.1)

                # get bisectors
                top = 0.9
                v_grid_cpu_bc, ccf_cpu_bc = GRASS.calc_ccf(chunck_vac_wav ./ ((-bc_ms)/GRASS.c_ms + 1), chunk_flux_full, lines, [maximum(chunk_flux_full) - minimum(chunk_flux_full)], resolution)
                vel_neid, int_neid = GRASS.calc_bisector(v_grid_cpu_bc, ccf_cpu_bc, nflux=50, top=top)
                ax2.plot(vel_neid[4:end], int_neid[4:end], c=colors[1])
                ax2.scatter(vel_neid[4:end], int_neid[4:end], c=colors[1])
                # set limits
                ax2.set_ylim(minimum(chunk_flux_full) - 0.05, 1.05)
                ax2.set_yticks([])

                idx = findfirst('_', line_names[i])
                fig.suptitle("$(line_names[i][1:idx-1]) $(round(vacwav[i]; digits = 1))", y=0.95)
                ax1.set_ylabel("Normalized Flux")
                ax1.set_xlabel("Wavelength (Å)")
                ax2.set_xlabel("Radial Velocity (m/s)")
                fig.savefig("eclipse_figures/Spectrum/timestamp_$(j).png")
            end
        end
        RV_all_lines[i] = RV_list
        RV_error_all_lines[i] = RV_error_list
    end
    return RV_all_lines, RV_error_all_lines
end

# october
# RV_all_lines, RV_error_all_lines = line_rvs_ccf(line_names, vacwav, orders, timestamps_october, path_october)
# @save "/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2"
# jldopen("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/neid_RVlinebyline.jld2", "a+") do file
#     file["name"] = line_names 
#     file["rv"] = RV_all_lines 
#     file["rv_error"] = RV_error_all_lines 
# end

# #for last timestamp (out of transit) of October eclipse, line comparsion between IAG, GRASS, and NEID
# last_timestamp_lines(line_names, airwav, vacwav, orders, neid_timestamps_october, timestamps_october, path_october, "neid_october_N_50", "KSSD")

#for last timestamp (out of transit) of October eclipse, line and bisector comparsion between GRASS and NEID
GRASS_comparison(line_names, airwav, vacwav, orders, neid_timestamps_october, timestamps_october, path_october, "neid_october_N_50", "KSSD")