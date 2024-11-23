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

# determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)

# line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]

datadir = string(abspath("data"))
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

# NEID location
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

function line_comp(line_names, airwav, vacwav, orders, neid_timestamps, obs_lat, obs_long, alt,
                                timestamps, path, filename, LD_type)
    # convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    variable_names = ["zenith", "dA_total_proj", "idx1", "idx3", "mu_grid", "z_rot_sub", "mu", "ax_codes", "dA", "N"]

    # Open the JLD2 file and read the variables into a dictionary
    data = jldopen("data/solar_disk/$(filename).jld2", "r") do file
        Dict(var => read(file, var) for var in variable_names)
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

    N = data["N"]
    Nt = length(time_stamps)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

    full_pixels = 2048:7168
    for i in 1:length(line_names)
        order_index = orders[i]
        min_wav = vacwav[i] - 2
        max_wav = vacwav[i] + 2

        # get the depth for the simulation
        sim_depth = df[i, "optimized_depth"]
        # simulate the spectrum
        lines = [airwav[i]]
        depths = [sim_depth]
        templates = [line_names[i]]
        resolution = 7e5
        spec = GRASS.SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution, oversampling=4.0)
    
        # simulate the spectrum 
        wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_type, obs_long, obs_lat, alt, 
                                                                time_stamps, zenith_mean, dA_total_proj, idx1, idx3, 
                                                                mu_grid, z_rot_sub, stored_μs, stored_ax_codes, stored_dA, 1.0, 
                                                                ext_toggle = false, verbose=true, use_gpu=false)

        f = FITS(joinpath(path, timestamps))
        header = read_header(f[1])
        bc_ms = (header["SSBRV0$order_index"]) * 1000                                                        

        # convolve GRASS spectrum to NEID resolution
        wavs_sim, flux_sim = GRASS.convolve_gauss(wavs_sim, flux_sim, new_res=11e4, oversampling=4.0)
        wavs_sim .= λ_air_to_vac.(wavs_sim)
        wavs_sim = wavs_sim ./ ((-636)/GRASS.c_ms + 1)

        # IAG flux 
        wavs_iag0, flux_iag0 = read_iag_atlas(isolate=true, airwav=airwav[i], buffer=2.0)
        flux_iag0 .-= 1 - maximum(flux_sim)

        spectrum = NEID.read_data(joinpath(path, timestamps); normalization = :blaze)
        # find where line is - NEID flux
        chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
        pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
        chunk_flux_full_o = chunk.flux[pixels]
        chunk_flux_full = chunk_flux_full_o ./ maximum(chunk_flux_full_o)
        chunk_flux_full .-= 1 - maximum(flux_sim)
        chunck_vac_wav = chunk.λ[pixels] ./ ((-bc_ms-636)/GRASS.c_ms + 1)

        # correct vacwav for BC (already corrected for grav. redshit)
        vacwav[i] = vacwav[i] ./ ((-bc_ms)/GRASS.c_ms + 1)

        # plot spectrum
        plt.figure()
        plt.plot(wavs_iag0, flux_iag0, label = "IAG", color = "g")
        plt.plot(chunck_vac_wav, chunk_flux_full, label = "NEID", color = "b")
        plt.scatter(chunck_vac_wav, chunk_flux_full, color = "b", s = 5)
        plt.plot(wavs_sim, flux_sim, label = "GRASS", color = "r")
        plt.scatter(wavs_sim, flux_sim, color = "r", s = 5)
        plt.axvline(x = vacwav[i])
        plt.legend()

        pixels = NEID.find_pixels_for_line_in_chunk(chunk, vacwav[i] - 0.25, vacwav[i] + 0.25)
        chunk_flux = chunk.flux[pixels]
        chunk_flux ./= maximum(chunk_flux_full_o)
        chunk_flux .-= 1 - maximum(flux_sim)
        chunck_vac_wav = chunk.λ[pixels] ./ ((-bc_ms-636)/GRASS.c_ms + 1)

        sim_shift_wl = (round(wavs_sim[findmin(flux_sim)[2]] - vacwav[i]; digits = 6))
        data_shift_wl = (round(chunck_vac_wav[findmin(chunk_flux)[2]] - vacwav[i]; digits = 6))
        sim_shift_vel = round(sim_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)
        data_shift_vel = round(data_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)

        plt.text(wavs_sim[120], 0.55, "Sim. Shift $sim_shift_wl Å")
        plt.text(wavs_sim[120], 0.5, "Sim. Shift $sim_shift_vel m/s")
        plt.text(wavs_sim[120], 0.45, "Data Shift $data_shift_wl Å")
        plt.text(wavs_sim[120], 0.4, "Data Shift $data_shift_vel m/s")

        plt.title("air wav: $(line_names[i])")
        plt.gca().invert_xaxis()
        plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/tests/NEID_spectrum/Spectrum_comp/spectra_$(i).png")
    end
end

function GRASS_comparison(line_names, airwav, vacwav, orders, neid_timestamps, obs_lat, obs_long, alt,
                            timestamps, path, filename, LD_type)
    # convert from utc to et as needed by SPICE
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
    # iterate through lines and determine line RV for eclipse EM curve
    for i in 1:length(line_names)

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

        # get the depth for the simulation
        sim_depth = df[i, "optimized_depth"]
        # simulate the spectrum
        lines = [airwav[i]]
        depths = [sim_depth]
        templates = [line_names[i]]
        resolution = 7e5
        spec = SpecParams(lines=lines, depths=depths, templates=templates, resolution=resolution, oversampling=4.0)
    
        # simulate the spectrum 
        wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, lines, LD_type, obs_long, obs_lat, alt, 
                                                                time_stamps, zenith_mean, dA_total_proj, idx1, idx3, 
                                                                mu_grid, z_rot_sub, stored_μs, stored_ax_codes, stored_dA, 1.0, 
                                                                ext_toggle = false, verbose=true, use_gpu=false)
        # convolve GRASS spectrum to NEID resolution
        wavs_sim, flux_sim = GRASS.convolve_gauss(wavs_sim, flux_sim, new_res=11e4, oversampling=4.0)
        wavs_sim .= λ_air_to_vac.(wavs_sim)
        wavs_sim = wavs_sim ./ ((-636)/GRASS.c_ms + 1)

        # convert from airwav to vacuum wavelength 
        spec.lambdas .= λ_air_to_vac.(spec.lambdas)
        spec.lines .= λ_air_to_vac.(spec.lines)

        f = FITS(joinpath(path, timestamps))
        header = read_header(f[1])
        bc_ms = (header["SSBRV0$order_index"]) * 1000

        # correct vacwav for BC (already corrected for grav. redshit)
        vacwav[i] = vacwav[i] ./ ((-bc_ms)/GRASS.c_ms + 1)

        spectrum = NEID.read_data(joinpath(path, timestamps); normalization = :blaze)
        # find where line is - NEID flux
        chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
        pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
        chunk_flux_full_o = chunk.flux[pixels]
        chunk_flux_full = chunk_flux_full_o ./ maximum(chunk_flux_full_o)
        chunk_flux_full .-= 1 - maximum(flux_sim)
        chunck_vac_wav = chunk.λ[pixels] ./ ((-bc_ms-636)/GRASS.c_ms + 1)

        # interpolate neid on synth wavelength grid
        itp = GRASS.linear_interp(chunck_vac_wav, chunk_flux_full, bc=NaN)
        chunk_flux_full = itp.(wavs_sim)
        chunck_vac_wav = copy(wavs_sim)

        v_grid_cpu_neid, ccf_cpu_neid = GRASS.calc_ccf(chunck_vac_wav, chunk_flux_full, spec)
        v_grid_cpu_sim, ccf_cpu_sim = GRASS.calc_ccf(wavs_sim, flux_sim, spec)

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
                ax1.plot(wavs_sim, flux_sim, marker="o", c="black", ms=3.0, lw=1.5, markevery=2, label="Synthetic")
                ax1.plot(chunck_vac_wav, chunk_flux_full, marker="s", c=colors[1], ms=2.0, lw=1.0, markevery=2, label="neid")
                # plot resids
                ax2.plot(wavs_sim, chunk_flux_full .- flux_sim, c=colors[1], marker="s", ms=2.0, lw=0, markevery=2)
            
                # find limits
                idx_min = argmin(flux_sim)[1]
                idx1 = idx_min - findfirst(x -> x .> 0.95, flux_sim[idx_min:-1:1])
                idx2 = idx_min + findfirst(x -> x .> 0.95, flux_sim[idx_min:end])
            
                # set limits
                min_idx = argmin(chunk_flux_full)
                ax1.set_xlim(wavs_sim[idx1-50], wavs_sim[idx2+50])
                ax1.set_ylim(minimum(flux_sim) - 0.1, 1.1)
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
                sim_shift_wl = (round(wavs_sim[findmin(flux_sim)[2]] - vacwav[i]; digits = 6))
                data_shift_wl = (round(chunck_vac_wav[idx_min-10:idx_min+10][findmin(chunk_flux_full[idx_min-10:idx_min+10])[2]] - vacwav[i]; digits = 6))
                sim_shift_vel = round(sim_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)
                data_shift_vel = round(data_shift_wl * GRASS.c_ms / vacwav[i]; digits = 2)
    
                ax1.text(wavs_sim[120], 0.75, "Sim. Shift $sim_shift_wl Å")
                ax1.text(wavs_sim[120], 0.70, "Sim. Shift $sim_shift_vel m/s")
                ax1.text(wavs_sim[120], 0.65, "Data Shift $data_shift_wl Å")
                ax1.text(wavs_sim[120], 0.60, "Data Shift $data_shift_vel m/s")
            
                # save the plot
                fig.subplots_adjust(wspace=0.05)
                fig.savefig(joinpath("/storage/home/efg5335/work/Eclipse_GRASS/tests/NEID_spectrum/GRASS_comp", line_name * "_line.png"))
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
                ax1.set_ylim(minimum(flux_sim) - 0.05, 1.05)
                ax2.set_xlim(-35, 35)
                ax2.set_ylim(minimum(flux_sim) - 0.05, 1.05)
                ax1.set_xlabel(L"{\rm Relative\ Velocity\ (m\ s^{-1})}", fontsize=13)
                ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
                ax2.set_xlabel(L"{\rm neid\ -\ GRASS\ (m\ s^{-1})}", fontsize=13)
                ax1.legend(labelspacing=0.25)
                ax2.legend(loc="upper right", prop=Dict("size"=>12.5), labelspacing=0.25)
            
                # set the title
                fig.suptitle(("\${\\rm " * title * "}\$"), y=0.95)
            
                # save the plot
                fig.subplots_adjust(hspace=0.05)
                fig.savefig(joinpath("/storage/home/efg5335/work/Eclipse_GRASS/tests/NEID_spectrum/GRASS_comp", line_name * "_bisector.png"))
                plt.clf(); plt.close()
                return nothing
            end
        comparison_plots()
    end
end

function line_rvs_ccf(line_names, vacwav, orders, timestamps, path)
    full_pixels = 2048:7168

    resolution = 11e4

    RV_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)
    RV_error_all_lines = Vector{Vector{Float64}}(undef,length(line_names)...)
    # iterate through lines and determine line RV for eclipse EM curve
    Threads.@threads for i in 1:length(line_names)
        order_index = orders[i]
        min_wav = vacwav[i] - 2
        max_wav = vacwav[i] + 2

        lines = [vacwav[i]]

        RV_list = Vector{Float64}(undef,length(timestamps)...)
        RV_error_list = Vector{Float64}(undef,length(timestamps)...)
        Threads.@threads for j in 1:length(timestamps)
            f = FITS(joinpath(path, timestamps[j]))
            # header = read_header(f[1])
            # bc_ms = (header["SSBRV0$order_index"]) * 1000

            spectrum = NEID.read_data(joinpath(path, timestamps[j]); normalization = :blaze)
            # find where line is - NEID flux
            chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
            pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
            chunk_flux_full = chunk.flux[pixels]
            chunk_flux_full ./= maximum(chunk_flux_full)
            chunck_vac_wav = chunk.λ[pixels]
            chunck_var = chunk.var[pixels] ./ maximum(chunk_flux_full)^2

            v_grid_cpu, ccf_cpu, ccf_var_out = GRASS.calc_ccf(chunck_vac_wav, chunk_flux_full, chunck_var, lines, [maximum(chunk_flux_full) - minimum(chunk_flux_full)], resolution)
            rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu, ccf_var_out)  

            RV_list[j] = rvs_cpu
            RV_error_list[j] = sigs_cpu

            # if line_names[i] == "FeI_5383"
            #     fig = plt.figure(figsize=(6.4,4.8))
            #     gs = mpl.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1.1], figure=fig, wspace=0.05)
            #     ax1 = fig.add_subplot(gs[1])
            #     ax2 = fig.add_subplot(gs[2])

            #     ax1.plot(chunck_vac_wav ./ ((-bc_ms)/GRASS.c_ms + 1), chunk_flux_full, c=colors[1])
            #     ax1.scatter(chunck_vac_wav ./ ((-bc_ms)/GRASS.c_ms + 1), chunk_flux_full, c=colors[1])
            #     ax1.axvline(x = vacwav[i], c=colors[1])
            #     # find limits
            #     idx_min = argmin(chunk_flux_full)
            #     idx1 = idx_min - findfirst(x -> x .> 0.95, chunk_flux_full[idx_min:-1:1])
            #     idx2 = idx_min + findfirst(x -> x .> 0.95, chunk_flux_full[idx_min:end])
            #     # set limits
            #     ax1.set_xlim(chunck_vac_wav[idx1-20], chunck_vac_wav[idx2+20])
            #     ax1.set_ylim(minimum(chunk_flux_full) - 0.1, 1.1)

            #     # get bisectors
            #     top = 0.9
            #     v_grid_cpu_bc, ccf_cpu_bc = GRASS.calc_ccf(chunck_vac_wav ./ ((-bc_ms)/GRASS.c_ms + 1), chunk_flux_full, lines, [maximum(chunk_flux_full) - minimum(chunk_flux_full)], resolution)
            #     vel_neid, int_neid = GRASS.calc_bisector(v_grid_cpu_bc, ccf_cpu_bc, nflux=50, top=top)
            #     ax2.plot(vel_neid[4:end], int_neid[4:end], c=colors[1])
            #     ax2.scatter(vel_neid[4:end], int_neid[4:end], c=colors[1])
            #     # set limits
            #     ax2.set_ylim(minimum(chunk_flux_full) - 0.05, 1.05)
            #     ax2.set_yticks([])

            #     idx = findfirst('_', line_names[i])
            #     fig.suptitle("$(line_names[i][1:idx-1]) $(round(vacwav[i]; digits = 1))", y=0.95)
            #     ax1.set_ylabel("Normalized Flux")
            #     ax1.set_xlabel("Wavelength (Å)")
            #     ax2.set_xlabel("Radial Velocity (m/s)")
            #     fig.savefig("/storage/home/efg5335/work/Eclipse_GRASS/tests/Spectrum/Fe_5384_movie/timestamp_$(j)_new.png")
            # end
        end
        RV_all_lines[i] = RV_list
        RV_error_all_lines[i] = RV_error_list
    end
    return RV_all_lines, RV_error_all_lines
end

function neid_eclipse()
    path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
    timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
    timestamps_october = timestamps_full_october[16:length(timestamps_full_october)-150]

    RV_all_lines, RV_error_all_lines = line_rvs_ccf(line_names, vacwav, orders, timestamps_october, path_october)
    @save "neid_RVlinebyline.jld2"
    jldopen("neid_RVlinebyline.jld2", "a+") do file
        file["name"] = line_names 
        file["rv"] = RV_all_lines 
        file["rv_error"] = RV_error_all_lines 
    end
end

function neid_nxt_day()
    path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/15/"
    files_and_dirs = readdir(path_october)
    files = filter(f -> isfile(joinpath(path_october, f)), files_and_dirs)
    timestamps_october = files[4:length(files)-145]

    RV_all_lines, RV_error_all_lines = line_rvs_ccf(line_names, vacwav, orders, timestamps_october, path_october)
    @save "neid_RVlinebyline_nxt_day.jld2"
    jldopen("neid_RVlinebyline_nxt_day.jld2", "a+") do file
        file["name"] = line_names 
        file["rv"] = RV_all_lines 
        file["rv_error"] = RV_error_all_lines 
    end
end

function neid_nxt_day()
    path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/15/"
    files_and_dirs = readdir(path_october)
    files = filter(f -> isfile(joinpath(path_october, f)), files_and_dirs)
    timestamps_october = files[4:length(files)-145]

    RV_all_lines, RV_error_all_lines = line_rvs_ccf(line_names, vacwav, orders, timestamps_october, path_october)
    @save "neid_RVlinebyline_nxt_day.jld2"
    jldopen("neid_RVlinebyline_nxt_day.jld2", "a+") do file
        file["name"] = line_names 
        file["rv"] = RV_all_lines 
        file["rv_error"] = RV_error_all_lines 
    end
end

function compare_spec_sources() # out of transit line comparsion between IAS, GRASS, and NEID
    neid_timestamps = ["2023-10-14T19:03:06.500000"] # last timestamp (out of transit)
    path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
    timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
    timestamps = timestamps_full_october[length(timestamps_full_october)-150]

    line_comp(line_names, airwav, vacwav, orders, neid_timestamps, obs_lat, obs_long, alt,
                timestamps, path_october, "neid_october_N_50", "SSD")

    GRASS_comparison(line_names, airwav, vacwav, orders, neid_timestamps, obs_lat, obs_long, alt,
                            timestamps, path_october, "neid_october_N_50", "SSD")
end

# neid_eclipse()
# neid_nxt_day()
compare_spec_sources()