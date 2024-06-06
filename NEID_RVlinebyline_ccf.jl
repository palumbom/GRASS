# 1. air wavelength to vaccuum
# 2. grass rv in neid resolution / pixel count
#     make line plot overploting grass and neid lines as scatters
# 3. calculate sun and WD g. redshift
# 4. individidual line in mask rather than spectra
# 5. get timing right for line by line
# 6. get RV error (already given) should match residual scatter rms out fo eclipse
# 7. turn on convective blueshift, turn off variability instead
# 8. once above complete, run one line at full resolution

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

refractive_index = 1.0003

#determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
#line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
min_wav_list = [1.3, 1.3, 1, 1.3, 2, 1.1, 1.2, 1.3, 1.1, 1.3, 1.3, 1.35, 1.1, 1.2, 1.35, 1.35, 1.35, 1.4, 1.35, 1.35, 1.35, 1.35]
max_wav_list = [1.6, 1.7, 2, 2, 2, 2, 1.7, 1.7, 1.9, 1.7, 1.7, 1.7, 2, 2, 1.9, 1.9, 2, 2, 1.9, 1.9, 2.1, 1.95]

#determine timestamps
#October Eclipse
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
#iterate through lines and determine line RV for eclipse EM curve
for i in 1:length(line_names)
    if i == 5
        continue
    end

    order_index = orders[i]
    min_wav = airwav[i] - 2
    max_wav = airwav[i] + 2

    #get the depth for the simulation
    sim_depth = df[i, "optimized_depth"]
    #simulate the spectrum
    lines = [airwav[i]]
    depths = [sim_depth]
    templates = [line_names[i]]
    resolution = 7e5
    spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution, buffer=1.5, oversampling=2.0)

    RV_list = Vector{Float64}(undef,length(timestamps)...)
    for j in 1:length(timestamps)

        f = FITS(joinpath(path, timestamps[j]))
        header = read_header(f[1])
        bc_z = (header["SSBZ0$order_index"])

        spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)

        #find where line is
        chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
        pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav) 
        chunck_air_wav = (chunk.λ[pixels] + bc_z * chunk.λ[pixels]) / refractive_index
        chunk_flux = chunk.flux[pixels]
        continuum = maximum(chunk.flux[pixels])

        # #plot specific line of interest
        # plt.figure()
        # plt.plot(chunck_air_wav, chunk_flux)

        pixels = NEID.find_pixels_for_line_in_chunk(chunk, airwav[i] + min_wav_list[i], airwav[i] + max_wav_list[i]) 
        chunck_air_wav = (chunk.λ[pixels] + bc_z * chunk.λ[pixels]) / refractive_index
        chunk_flux = chunk.flux[pixels]

        min_flux_ind = argmin(chunk_flux)
        itp1 = GRASS.linear_interp(chunck_air_wav, chunk_flux)
        wav_range = LinRange(chunck_air_wav[min_flux_ind] - 0.5, chunck_air_wav[min_flux_ind] + 0.5, 100)
        
        # plt.plot(wav_range, itp1.(wav_range), color = "purple")
        # plt.plot(chunck_air_wav, chunk_flux)
        # plt.axvline(x=airwav[i])
        # plt.title(line_names[i])
        # plt.gca().invert_xaxis()
        # plt.savefig("spectra_zoom_$i.png")

        # v_grid_cpu, ccf_cpu = GRASS.calc_ccf(wav_range, itp1.(wav_range), spec)
        # rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

        # RV_list[j] = rvs_cpu
    end

    # RV_all_lines[i] = RV_list

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
# end