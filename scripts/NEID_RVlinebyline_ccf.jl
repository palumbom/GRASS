using CSV
using GRASS
using SPICE
using FITSIO
using PyPlot
using FileIO
using JLD2
using Statistics
using DataFrames
using EchelleCCFs
using RvSpectMLBase
using EchelleInstruments

neid_timestamps = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
#convert from utc to et as needed by SPICE
time_stamps = utc2et.(neid_timestamps)
#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

refractive_index = 1.0003

#determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
airvac = airwav .* refractive_index
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

datadir = string(abspath("data"))
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

full_pixels = 2048:7168
RV_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)
RV_error_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)
#iterate through lines and determine line RV for eclipse EM curve
for i in 1:length(line_names)

    order_index = orders[i]
    min_wav = airvac[i] - 2
    max_wav = airvac[i] + 2

    #get the depth for the simulation
    sim_depth = df[i, "optimized_depth"]
    #simulate the spectrum
    lines = [airwav[i]]
    depths = [sim_depth]
    templates = [line_names[i]]
    resolution = 11e4
    spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution, oversampling=4.0) #buffer=1.5,
    
    N = 50
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    # simulate the spectrum
    wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, obs_long, obs_lat, alt, lines ./ 10.0, time_stamps, "Optical", verbose=true, use_gpu=false)

    #convert from airwav to vacuum wavelength 
    spec.lambdas .= spec.lambdas * refractive_index
    spec.lines .= spec.lines * refractive_index
    
    RV_list = Vector{Float64}(undef,length(timestamps)...)
    RV_error_list = Vector{Float64}(undef,length(timestamps)...)
    for j in 1:1#length(timestamps)
        spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)

        #find where line is
        chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
        pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
        chunk_flux = chunk.flux[pixels]
        chunk_flux ./= maximum(chunk_flux)

        #plot spectrum
        plt.figure()
        plt.plot(chunk.λ[pixels], chunk_flux, color = "b")
        plt.scatter(chunk.λ[pixels], chunk_flux, label = "NEID", color = "b", s = 5)
        plt.plot(wavs_sim, flux_sim[((j-1)*length(wavs_sim) + 1):(length(wavs_sim)*j)], color = "r")
        plt.scatter(wavs_sim, flux_sim[((j-1)*length(wavs_sim)+ 1):(length(wavs_sim)*j)], label = "GRASS", color = "r", s = 5)
        plt.legend()
        plt.axvline(x=airvac[i])
        plt.title(line_names[i])
        plt.gca().invert_xaxis()
        plt.savefig("spectra_zoom_$i.png")

        # v_grid_cpu, ccf_cpu = GRASS.calc_ccf(chunk.λ[pixels], chunk_flux, spec)
        # rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        # RV_list[j] = rvs_cpu
        # RV_error_list[j] = sigs_cpu
    end

    # RV_all_lines[i] = RV_list
    # RV_error_all_lines[i] = RV_error_list

    # #plot rm curve for specific lines
    # plt.figure()
    # plt.plot(1:length(timestamps), RV_all_lines[i])
    # plt.title(airwav[i])
    # plt.savefig("rm_$i.png")
end

# @save "neid_RVlinebyline.jld2"
# jldopen("neid_RVlinebyline.jld2", "a+") do file
#     file["name"] = line_names 
#     file["rv"] = RV_all_lines 
#     file["rv_error"] = RV_error_all_lines 
# end