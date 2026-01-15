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
using Interpolations
using EchelleCCFs: λ_air_to_vac, λ_vac_to_air
using Dierckx
# plotting
using LaTeXStrings
import PyPlot
plt = PyPlot
mpl = plt.matplotlib
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

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

datadir = string(abspath("data"))
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

# determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)
λrest = GRASS.get_rest_wavelength(lp)
depth = GRASS.get_depth(lp)
lfile = GRASS.get_file(lp)

# line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]

variable_names = ["lambda_min", "lambda_max", "flux_min", "pixel_mean", "pixels", "neid_wavelength"]
ext_coeff_array = [0.15452995224327976, 0.15256098077094832, 0.14720055859068512, 0.154895798933504, 0.15181381895180662, 0.15107508233588227, 0.15116772762156633, 0.14882114581650618, 0.14865707189399568, 0.1494903120065096, 0.16011027092744037, 0.15593033972594958, 0.14195968590211427, 0.15401904166429853, 0.1277699772941639, 0.12709315507233226, 0.12820346527304866, 0.11702310600015708, 0.1435320747844216, 0.12380490304619193, 0.12450734135297492, 0.12101777355247835]

# NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

GRASS.get_kernels()

function resolving_power(wavelengths, intensities)
    interpolation = CubicSpline(wavelengths, intensities)

    x_fine = wavelengths[1]:0.01:wavelengths[length(wavelengths)]
    y_fine = interpolation.(x_fine)
    y_fine = [x[] for x in y_fine]

    peak_index = argmin(y_fine)
    peak_wavelength = x_fine[peak_index]
    
    min_intensity = minimum(y_fine)
    I_half = (1 + min_intensity) / 2
    dists = abs.(y_fine .- I_half)
    # Get indices of the two closest values
    closest_indices = sortperm(dists)[1:2]
    fwhm = abs((x_fine[closest_indices[1]]) - (x_fine[closest_indices[2]]))

    resolving_power = peak_wavelength / fwhm 
    return resolving_power, peak_wavelength
end

function line_rv(line_names, airwav, vacwav, orders, neid_timestamps, obs_lat, obs_long, alt, timestamps, path)
    RV_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)
    RV_error_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)

    RV_list_neid = Vector{Float64}(undef,length(timestamps)...)
    RV_err_list_neid = Vector{Float64}(undef,length(timestamps)...)

    lambda_min = Vector{Float64}(undef,length(timestamps)...)
    lambda_max = Vector{Float64}(undef,length(timestamps)...)
    pixel_mean = Vector{Float64}(undef,length(timestamps)...)
    flux_min = Vector{Float64}(undef,length(timestamps)...)
    pixels = Vector{Vector{Float64}}(undef,length(timestamps)...)
    neid_wavelength = Vector{Vector{Float64}}(undef,length(timestamps)...)

    full_pixels = 2048:7168 
    for i in 1:length(line_names)
        RV_list_neid = Vector{Float64}(undef,length(timestamps)...)
        RV_err_list_neid = Vector{Float64}(undef,length(timestamps)...)

        if i in vcat(10:22)
            Δv_max=6000.0
        elseif i in vcat(3:9)
            Δv_max=7000.0
        else
            Δv_max=8000.0
        end

        for j in 1:length(timestamps)
            order = orders[i]
            min_wav = vacwav[i] - 5.0
            max_wav = vacwav[i] + 5.0
            buffer = 0.2

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            chunk = NEID.ChunkOfSpectrum(spectrum, order, full_pixels) 
            pixel_normal = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            pixels_inner = NEID.find_pixels_for_line_in_chunk(chunk, vacwav[i] - buffer, vacwav[i] + buffer)
            pixel = mean(pixels_inner)
            order_wavelengths = chunk.λ[pixels_inner]
            chunk_flux = chunk.flux[pixels_inner]
            chunk_flux = chunk_flux ./ maximum(chunk.flux[pixel_normal])

            # f = FITS(joinpath(path, timestamps[j]))
            # hdr = read(f[11])
            # wav = read(f[8])[:, order]
            # indices_in_range = findall(x -> min_wav <= x <= max_wav, wav)
            # spl = Spline1D(wav[indices_in_range], hdr[1, indices_in_range, order]; k=3)
            # chunk_flux = chunk_flux ./ spl(order_wavelengths)

            v_grid, ccf = GRASS.calc_ccf(order_wavelengths, chunk_flux, [vacwav[i]], [1.0 - minimum(chunk_flux)], 11e4, Δv_max=Δv_max)
            rvs, sigs = GRASS.calc_rvs_from_ccf(v_grid, ccf)  
            RV_list_neid[j] = rvs
            RV_err_list_neid[j] = sigs

            lambda_min[j] = order_wavelengths[1]                                                 
            lambda_max[j] = order_wavelengths[length(order_wavelengths)] 
            pixels[j] = Array(pixels_inner)
            neid_wavelength[j] = order_wavelengths
            pixel_mean[j] = pixel
            flux_min[j] = minimum(chunk_flux)
        end

        RV_all_lines[i] = RV_list_neid
        RV_error_all_lines[i] = RV_err_list_neid

        @save "data/NEID_convolution_info$(i).jld2"
        jldopen("data/NEID_convolution_info$(i).jld2", "a+") do file
            file["lambda_min"] = deepcopy(lambda_min)
            file["lambda_max"] = deepcopy(lambda_max) 
            file["flux_min"] = deepcopy(flux_min) 
            file["pixel_mean"] = deepcopy(pixel_mean)
            file["pixels"] = deepcopy(pixels) 
            file["neid_wavelength"] = deepcopy(neid_wavelength)
        end
    end
return RV_all_lines, RV_error_all_lines
end 

function line_comp(line_names, airwav, vacwav, orders, neid_timestamps, obs_lat, obs_long, alt, timestamps, path, LD_type)
    # convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    N = 197
    Nt = length(neid_timestamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=40)
    resolution = 7e5

    full_pixels = 2048:7168 
    for i in 1:length(line_names)
        if i in vcat(10:22)
            Δv_max=6000.0
        elseif i in vcat(3:9)
            Δv_max=7000.0
        else
            Δv_max=8000.0
        end

        println("\t>>> Template: " * string(splitdir(lfile[i])[2]))
        data = jldopen("data/NEID_convolution_info$(i).jld2", "r") do file
            Dict(var => read(file, var) for var in variable_names)
        end
        lambda_min = data["lambda_min"]
        lambda_max = data["lambda_max"]
        flux_min = data["flux_min"]
        pixel_mean = data["pixel_mean"]
        pixels = data["pixels"]
        neid_wavelength = data["neid_wavelength"]

        # set up parameters for synthesis
        lines = [λrest[i]]
        templates = [lfile[i]]
        depths = [df[i, "optimized_depth"]]
        variability = trues(length(lines))  # whether or not the bisectors should "dance"
        blueshifts = zeros(length(lines))   # set convective blueshift value

        # make the spec composite type instances
        spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
                        blueshifts=blueshifts, templates=templates, resolution=resolution) 

        lambdas, outspec = GRASS.synth_Eclipse_gpu(spec, disk, true, Float64, falses(disk.Nt), LD_type, obs_long, obs_lat, alt, time_stamps, lines, ext_coeff_array[i], true, true)
        lambdas .= λ_air_to_vac.(lambdas)

        for j in 1:length(time_stamps)
            outspec_t = outspec[:, j]
            outspec_t = outspec_t ./ maximum(outspec_t)                                                            

            model_pixels, flux_sim, dlam_dpix = neid_asymmetric_lsf.lsf_model_export(orders[i], round(Int, pixel_mean[j]), lambdas, outspec_t, vacwav[i])
            model_wavelength = (dlam_dpix .* model_pixels) .+ vacwav[i]

            order = orders[i]
            min_wav = vacwav[i] - 5.0
            max_wav = vacwav[i] + 5.0
            buffer = 0.2

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            chunk = NEID.ChunkOfSpectrum(spectrum, order, full_pixels) 
            pixel_normal = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            pixels_inner = NEID.find_pixels_for_line_in_chunk(chunk, vacwav[i] - buffer, vacwav[i] + buffer)
            order_wavelengths = chunk.λ[pixels_inner]
            chunk_flux = chunk.flux[pixels_inner]
            chunk_flux = chunk_flux ./ maximum(chunk.flux[pixel_normal])

            plt.figure()
            plt.plot(order_wavelengths, chunk_flux, label = "NEID R~$(round(resolving_power(order_wavelengths, chunk_flux)[1]))", color = "b")
            plt.scatter(order_wavelengths, chunk_flux, color = "b", s = 5)
            plt.plot(model_wavelength, flux_sim, label = "GRASS Asymmetric Convolved R~$(round(resolving_power(model_wavelength, flux_sim)[1]))", color = "g")
            plt.scatter(model_wavelength, flux_sim, color = "g", s = 5)
            plt.plot(lambdas, outspec_t, label = "GRASS R~$(round(resolving_power(lambdas, outspec_t)[1]))", color = "r")
            plt.scatter(lambdas, outspec_t, color = "r", s = 5)
            plt.axvline(x = vacwav[i], color = "k")
            # plt.axvline(x = resolving_power(order_wavelengths, chunk_flux)[2], color = "b")
            # plt.axvline(x = resolving_power(model_wavelength, flux_sim)[2], color = "g")
            plt.legend(loc="lower left", fontsize=8)
            plt.xlabel("Wavelength Å")
            plt.xlim(order_wavelengths[1], order_wavelengths[length(order_wavelengths)])
            plt.ylabel("Normalized Flux")
            plt.title("air wav: $(line_names[i])")
            plt.gca().invert_xaxis()
            plt.savefig("line_profiles/line_$(i).png")
        end
    end
end 

function compute_linebyline_neid_rv() 
    neid_timestamps = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
    path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
    timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
    timestamps = timestamps_full_october[16:length(timestamps_full_october)-150]

    RV_all_lines, RV_error_all_lines = line_rv(line_names, airwav, vacwav, orders, neid_timestamps[1:130], obs_lat, obs_long, alt,
                timestamps[1:130], path_october)

    @save "neid_RVlinebyline.jld2"
    jldopen("neid_RVlinebyline.jld2", "a+") do file
        file["name"] = line_names 
        file["rv"] = RV_all_lines 
        file["rv_error"] = RV_error_all_lines 
    end
end

function spectral_lines() 
    neid_timestamps = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
    path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
    timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
    timestamps = timestamps_full_october[16:length(timestamps_full_october)-150]

    line_comp(line_names, airwav, vacwav, orders, neid_timestamps[130:130], obs_lat, obs_long, alt,
                timestamps[130:130], path_october, "SSD_4parameter")
end

compute_linebyline_neid_rv()