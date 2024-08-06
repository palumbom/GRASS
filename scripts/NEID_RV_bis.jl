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
# plotting
using LaTeXStrings
import PyPlot
plt = PyPlot
mpl = plt.matplotlib
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

GRASS.get_kernels()

#determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
files = GRASS.get_file(lp)
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)

#line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]

datadir = string(abspath("data"))
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)
plotdir = string(abspath(joinpath("neid_figures", "October_zoom_in_comp_aligned")))
if !isdir(plotdir)
    mkdir(plotdir)
end

# decide whether to use gpu
use_gpu = false

neid_timestamps_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
neid_timestamps_april =["2024-04-08T20:11:26", "2024-04-08T20:10:03", "2024-04-08T20:08:41", "2024-04-08T20:07:18", "2024-04-08T20:05:55", "2024-04-08T20:04:33", "2024-04-08T20:03:10", "2024-04-08T20:01:47", "2024-04-08T20:00:25", "2024-04-08T19:59:02", "2024-04-08T19:57:39", "2024-04-08T19:56:16", "2024-04-08T19:54:54", "2024-04-08T19:53:31", "2024-04-08T19:52:08", "2024-04-08T19:50:46", "2024-04-08T19:49:23", "2024-04-08T19:48:00", "2024-04-08T19:46:38", "2024-04-08T19:45:15", "2024-04-08T19:43:52", "2024-04-08T19:42:30", "2024-04-08T19:41:07", "2024-04-08T19:39:44", "2024-04-08T19:38:22", "2024-04-08T19:36:59", "2024-04-08T19:35:36", "2024-04-08T19:34:14", "2024-04-08T19:32:51", "2024-04-08T19:31:28", "2024-04-08T19:30:06", "2024-04-08T19:28:43", "2024-04-08T19:27:20", "2024-04-08T19:25:57", "2024-04-08T19:24:35", "2024-04-08T19:23:12", "2024-04-08T19:21:49", "2024-04-08T19:20:27", "2024-04-08T19:19:04", "2024-04-08T19:17:41", "2024-04-08T19:16:19", "2024-04-08T19:14:56", "2024-04-08T19:13:33", "2024-04-08T19:12:11", "2024-04-08T19:10:48", "2024-04-08T19:09:25", "2024-04-08T19:08:03", "2024-04-08T19:06:40", "2024-04-08T19:05:17", "2024-04-08T19:03:55", "2024-04-08T19:02:32", "2024-04-08T19:01:09", "2024-04-08T18:59:47", "2024-04-08T18:58:24", "2024-04-08T18:57:01", "2024-04-08T18:55:38", "2024-04-08T18:54:16", "2024-04-08T18:52:53", "2024-04-08T18:51:30", "2024-04-08T18:50:08", "2024-04-08T18:48:45", "2024-04-08T18:47:22", "2024-04-08T18:46:00", "2024-04-08T18:44:37", "2024-04-08T18:43:14", "2024-04-08T18:41:52", "2024-04-08T18:40:29", "2024-04-08T18:39:06", "2024-04-08T18:37:44", "2024-04-08T18:36:21", "2024-04-08T18:34:58", "2024-04-08T18:33:36", "2024-04-08T18:32:13", "2024-04-08T18:30:50", "2024-04-08T18:29:28", "2024-04-08T18:28:05", "2024-04-08T18:26:42", "2024-04-08T18:25:19", "2024-04-08T18:23:57", "2024-04-08T18:22:34", "2024-04-08T18:21:11", "2024-04-08T18:19:49", "2024-04-08T18:18:26", "2024-04-08T18:17:03", "2024-04-08T18:15:41", "2024-04-08T18:14:18", "2024-04-08T18:12:55", "2024-04-08T18:11:33", "2024-04-08T18:10:10", "2024-04-08T18:08:47", "2024-04-08T18:07:25", "2024-04-08T18:06:02", "2024-04-08T18:04:39", "2024-04-08T18:03:17", "2024-04-08T18:01:54", "2024-04-08T18:00:31", "2024-04-08T17:59:09", "2024-04-08T17:57:46", "2024-04-08T17:56:23", "2024-04-08T17:55:01", "2024-04-08T17:53:38", "2024-04-08T17:52:15", "2024-04-08T17:50:52", "2024-04-08T17:49:30", "2024-04-08T17:48:07", "2024-04-08T17:46:44", "2024-04-08T17:45:22", "2024-04-08T17:43:59", "2024-04-08T17:42:36", "2024-04-08T17:41:14", "2024-04-08T17:39:51", "2024-04-08T17:38:28", "2024-04-08T17:37:06", "2024-04-08T17:35:43", "2024-04-08T17:34:20", "2024-04-08T17:32:58", "2024-04-08T17:31:35", "2024-04-08T17:30:12", "2024-04-08T17:28:50", "2024-04-08T17:27:27", "2024-04-08T17:26:04", "2024-04-08T17:24:42", "2024-04-08T17:23:19", "2024-04-08T17:21:56", "2024-04-08T17:20:33", "2024-04-08T17:19:11", "2024-04-08T17:17:48", "2024-04-08T17:16:25", "2024-04-08T17:15:03", "2024-04-08T17:13:40", "2024-04-08T17:12:17", "2024-04-08T17:10:55", "2024-04-08T17:09:32", "2024-04-08T17:08:09", "2024-04-08T17:06:47", "2024-04-08T17:05:24", "2024-04-08T17:04:01", "2024-04-08T17:02:39", "2024-04-08T17:01:16", "2024-04-08T16:59:53", "2024-04-08T16:58:31", "2024-04-08T16:57:08", "2024-04-08T16:55:45", "2024-04-08T16:54:23", "2024-04-08T16:53:00", "2024-04-08T16:51:37", "2024-04-08T16:50:14", "2024-04-08T16:48:52"]
# October Eclipse
path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
timestamps_october = timestamps_full_october[16:length(timestamps_full_october)-150]
#April Eclipse
path_april = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2024/04/08/"
timestamps_full_april = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_April/data/neid_april_data.csv", DataFrame)[!, "filename"]
timestamps_april = timestamps_full_april[101:length(timestamps_full_april)-10]

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938

function neid_bis(line_names, files, airwav, vacwav, orders, neid_timestamps, timestamps, path, filename, LD_type)
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    variable_names = ["zenith_mean", "dA_total_proj", "idx1", "idx3", "mu_grid", "z_rot_sub", "mu", "ax_codes", "dA", "N"]

    # Open the JLD2 file and read the variables into a dictionary
    data = jldopen("data/solar_disk/$(filename).jld2", "r") do file
        Dict(var => read(file, var) for var in variable_names)
    end

    N = data["N"]
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    full_pixels = 2048:7168
    # wavelength of line to synthesize/compare to neid
    for i in 1:length(files)
        if i == 11
            continue
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

        println(">>> Running " * line_names[i] * "...")

        # get properties from line
        line_name = line_names[i]
        order_index = orders[i]
        min_wav = vacwav[i] - 1.5
        max_wav = vacwav[i] + 1.5

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
                                                                dA_total_proj, idx1, idx3, mu_grid, z_rot_sub,
                                                                stored_μs, stored_ax_codes, stored_dA, [0.0], ext_toggle = false, verbose=true, use_gpu=false)
        # convolve GRASS spectrum to NEID resolution
        wavs_sim, flux_sim = GRASS.convolve_gauss(wavs_sim, flux_sim, new_res=11e4, oversampling=4.0)
        wavs_sim .= λ_air_to_vac.(wavs_sim)
        #convert from airwav to vacuum wavelength 
        spec.lambdas .= λ_air_to_vac.(spec.lambdas)
        spec.lines .= λ_air_to_vac.(spec.lines)

        for j in 1:1#length(timestamps)
            #GRASS flux for given timestamps
            flux_sim_j = flux_sim[((j-1)*length(wavs_sim) + 1):(length(wavs_sim)*j)]

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            #find where line is - NEID flux
            chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
            pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            wavs_neid = chunk.λ[pixels]
            flux_neid = chunk.flux[pixels]
            flux_neid ./= maximum(flux_neid)
    
            idxl = findfirst(x -> x .>= vacwav[i] - buff, wavs_neid)
            idxr = findfirst(x -> x .>= vacwav[i] + buff, wavs_neid)
            neid_bot = minimum(view(flux_neid, idxl:idxr))
            neid_depth = 1.0 - neid_bot
    
            # interpolate neid on synth wavelength grid
            itp = GRASS.linear_interp(wavs_neid, flux_neid, bc=NaN)
            flux_neid = itp.(wavs_sim)
            wavs_neid = copy(wavs_sim)

            v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim_j, spec)
            v_grid_neid, ccf_neid = GRASS.calc_ccf(wavs_neid, flux_neid, spec)

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
            v_grid_neid, ccf_neid = GRASS.calc_ccf(wavs_neid, flux_neid, spec)
            v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim_j, spec)

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
                ax1.plot(wavs_sim, flux_sim_j, marker="o", c="black", ms=3.0, lw=1.5, markevery=2, label="Synthetic")
                ax1.plot(wavs_neid, flux_neid, marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=2, label="neid")
                # plot resids
                ax2.plot(wavs_sim, flux_neid .- flux_sim_j, c=colors[1], marker="s", ms=2.0, lw=0, markevery=2)
            
                # find limits
                idx_min = argmin(flux_sim_j)
                idx1 = idx_min - findfirst(x -> x .> 0.95, flux_sim_j[idx_min:-1:1])
                idx2 = idx_min + findfirst(x -> x .> 0.95, flux_sim_j[idx_min:end])
            
                # set limits
                min_idx = argmin(flux_neid)
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
            
                # save the plot
                fig.subplots_adjust(wspace=0.05)
                fig.savefig(joinpath(plotdir, line_name * "_line.png"))
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
                fig.savefig(joinpath(plotdir, line_name * "_bisector.png"))
                plt.clf(); plt.close()
                return nothing
            end

            comparison_plots()
        end
    end
end

#for first timestamp (out of transit) of October eclipse, line and bisector comparsion between GRASS and NEID
neid_bis(line_names, files, airwav, vacwav, orders, neid_timestamps_october, timestamps_october, path_october, "neid_october_N_5", "KSSD")