using Pkg; Pkg.activate(".")
using CSV
using HDF5
using Glob
using LsqFit
using DataFrames
using Statistics
using Polynomials
using GRASS

import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# set LARS spectra absolute dir and read line info file
const data_dir = "/storage/group/ebf11/default/mlp95/lars_spectra/"
const line_info = CSV.read(GRASS.soldir * "line_info.csv", DataFrame)

# function to write line parameters file
function write_line_params(line_df::DataFrame; clobber::Bool=false)
    # get the filename
    line_dir = GRASS.soldir * line_df.name[1] * "/"
    prop_file =  line_dir * line_df.name[1] * "_line_properties.h5"

    # don't bother if the file already exists
    if isfile(prop_file) & !clobber
        println("\t >>> " * splitdir(prop_file)[end] * " already exists...")
        return nothing
    end

    # read in the IAG data and isolate the line
    iag_wavs, iag_flux = read_iag(isolate=true, airwav=line_df.air_wavelength[1])
    iag_depth = 1.0 - minimum(iag_flux)

    # write the line properties file
    println("\t >>> Writing " * line_df.name[1] * "_line_properties.h5")
    h5open(prop_file, "w") do fid
        create_group(fid, "properties")
        g = fid["properties"]

        # fill out attributes
        attr = attributes(g)
        for n in names(line_df)
            if ismissing(line_df[!, n][1])
                attr[n] = NaN
            else
                attr[n] = line_df[!, n][1]
            end
        end
        attr["depth"] = iag_depth
    end
    return nothing
end

function write_input_data(line_name, air_wavelength, fparams, wav, bis, dep, wid)
    # create output file name
    new_file = line_name * "_" * string(fparams[3]) * "_" * fparams[5] * "_" * fparams[6] * "_input.h5"

    # write the input data to the file
    h5open(GRASS.soldir * line_name * "/" * new_file, "w") do fid
        # create the group
        create_group(fid, "input_data")
        g = fid["input_data"]

        # fill out the datasets
        g["wavelengths"] = wav
        g["bisectors"] = bis
        g["depths"] = dep
        g["widths"] = wid

        # make attributes
        attr = attributes(g)
        attr["datetime"] = string(fparams[3])
        attr["air_wavelength"] = air_wavelength

        # convert mu to number
        mu_num = []
        for ch in fparams[5]
            push!(mu_num, tryparse(Int64, string(ch)))
        end

        new_string = prod(string.(mu_num[.!isnothing.(mu_num)]))
        if new_string[1] == '1'
            attr["mu"] = 1.0
            attr["axis"] = "c"
        elseif new_string[1] == '0'
            attr["mu"] = parse(Float64, "0." * new_string[2:end])
            attr["axis"] = fparams[6]
        end
    end
    return nothing
end

function fit_line_wings(wavs_iso, flux_iso)
    # get min and flux thresholds
    min = argmin(flux_iso)
    bot = flux_iso[min]
    dep = 1.0 - bot
    lidx50 = min - findfirst(x -> x .>= 0.5 * dep + bot, reverse(flux_iso[1:min]))
    ridx50 = findfirst(x -> x .>= 0.5 * dep + bot, flux_iso[min:end]) + min
    lidx60 = min - findfirst(x -> x .>= 0.6 * dep + bot, reverse(flux_iso[1:min]))
    ridx60 = findfirst(x -> x .>= 0.6 * dep + bot, flux_iso[min:end]) + min
    lidx70 = min - findfirst(x -> x .>= 0.7 * dep + bot, reverse(flux_iso[1:min]))
    ridx70 = findfirst(x -> x .>= 0.7 * dep + bot, flux_iso[min:end]) + min
    lidx80 = min - findfirst(x -> x .>= 0.8 * dep + bot, reverse(flux_iso[1:min]))
    ridx80 = findfirst(x -> x .>= 0.8 * dep + bot, flux_iso[min:end]) + min
    lidx90 = min - findfirst(x -> x .>= 0.9 * dep + bot, reverse(flux_iso[1:min]))
    ridx90 = findfirst(x -> x .>= 0.9 * dep + bot, flux_iso[min:end]) + min

    # isolate the line wings
    Δbot = 1
    core = min-Δbot:min+Δbot
    lwing = lidx90:lidx50
    rwing = ridx50:ridx90
    wavs_fit = vcat(wavs_iso[lwing], wavs_iso[core], wavs_iso[rwing])
    flux_fit = vcat(flux_iso[lwing], flux_iso[core], flux_iso[rwing])

    # set boundary conditions and initial guess
    # GOOD FOR FeI 5434 + others
    lb = [0.0, wavs_iso[min], 0.0, 0.0]
    ub = [1.0, wavs_iso[min], 0.5, 0.5]
    p0 = [1.0 - dep, wavs_iso[min], 0.02, 0.01]
    # GOOD FOR FeI 5434 + others

    # perform the fit
    fit = curve_fit(GRASS.fit_voigt, wavs_fit, flux_fit, p0, lower=lb, upper=ub)
    @show fit.param
    return fit
end

function preprocess_line(line_name::String; verbose::Bool=true)
    # create subdirectory structure if it doesn't already exist
    if !isdir(GRASS.soldir * line_name)
        mkdir(GRASS.soldir * line_name)
    end

    # find row with line info and write the line_params file
    line_df = subset(line_info, :name => x -> x .== line_name)
    write_line_params(line_df)

    # find all the spectra files associated with this line
    fits_files = Glob.glob("*.fits", data_dir * line_df.spectra_dir[1] * "/")

    # read in the spectrum and bin into 15-second bins
    # for f in fits_files
    for f in fits_files[1:1]
        # print the filename
        if verbose println("\t >>> Processing " * splitdir(f)[end]) end

        # get spec parameters
        fparams = GRASS.extract_line_params(f)
        wavs, flux = GRASS.bin_spectrum(GRASS.read_spectrum(f)...)
        plt.plot(wavs[:,1], flux[:,1]./maximum(flux[:,1]), c="k")

        # normalize the spectra
        flux ./= maximum(flux, dims=1)

        # loop over epochs in data and measure input data
        wav = zeros(100, size(wavs,2))
        bis = zeros(100, size(wavs,2))
        wid = zeros(100, size(wavs,2))
        dep = zeros(100, size(wavs,2))
        for t in 1:size(wavs, 2)
            # DEBUGGING STUFF
            if t != 1
                continue
            end
            # DEBUGGING STUFF

            # refine the location of the minimum
            idx = findfirst(x -> x .>= line_df.air_wavelength[1], wavs[:,t])
            min = argmin(flux[idx-50:idx+50, t]) + idx - 50
            bot = flux[min, t]
            dep = 1.0 - bot

            # isolate the line
            idx1 = min - findfirst(x -> x .>= 0.95, reverse(flux[1:min, t]))
            idx2 = findfirst(x -> x .>= 0.95, flux[min:end, t]) + min
            idxl = min - findfirst(x -> x .>= 0.9 * dep + bot, reverse(flux[1:min, t]))
            idxr = findfirst(x -> x .>= 0.9 * dep + bot, flux[min:end, t]) + min
            wavs_iso = copy(view(wavs, idx1:idx2, t))
            flux_iso = copy(view(flux, idx1:idx2, t))

            plt.axvline(wavs[idxl, t])
            plt.axvline(wavs[idxr, t])


            # fit the line wings
            fit = fit_line_wings(wavs_iso, flux_iso)

            # get fit curve on original grid of wavelengths
            flux_new = GRASS.fit_voigt(view(wavs, :, t), fit.param)

            # replace wings with model
            idx_lw = argmin(flux_new[1:min])
            idx_rw = argmin(flux_new[min:end])
            flux[1:idxl, t] .= flux_new[1:idxl]
            flux[idxr:end, t] .= flux_new[idxr:end]

            # do a quick "normalization"
            flux[:,t] ./= maximum(flux[:,t])

            # oversample it to get better precision in wings
            itp = GRASS.linear_interp(wavs[:,t], flux[:,t])
            wavs_meas = range(first(wavs[:,t]), last(wavs[:,t]), step=wavs[min,t]/1e6)
            flux_meas = itp.(wavs_meas)

            # TODO REVIEW BISECTOR CODE
            # measure a bisector
            # wav[:,t], bis[:,t] = GRASS.measure_bisector_interpolate(wavs_meas, flux_meas)
            # dep[:,t], wid[:,t] = GRASS.measure_width_interpolate(wavs_meas, flux_meas)

            # # set any widths less than zero to 0
            # idx = findall(x -> x .< 0.0, view(wid, :, t))
            # wid[idx,t] .= 0.0

            # # polyfit the last three points and extrapolate to continuum
            # idx8 = findfirst(x -> x .>= 0.8, dep[:,t])
            # idx9 = findfirst(x -> x .>= 0.9, dep[:,t])
            # pfit = Polynomials.fit(dep[idx9-10:idx9,t], wid[idx9-10:idx9,t], 2)
            # wid[idx9:end, t] = pfit.(dep[idx9:end, t])

            # DEBUGGING STUFF
            plt.plot(wavs[:,t], flux[:,t], c="tab:orange")
            plt.plot(wavs[:,t], flux_new, c="tab:purple")
            # plt.plot(wavs_iso, flux_iso, c="tab:blue")
            # plt.plot(wav[:,t], bis[:,t], c="k")
            plt.show()

            # plt.plot(dep[:,t], wid[:,t]); plt.show()

            # plt.plot(dep[1:idx9,t], wid[1:idx9,t], c="k", label="original")
            # plt.plot(dep[idx9:end,t], wid[idx9:end,t], label="extrap")
            # DEBUGGING STUFF
        end

        # DEBUGGING STUFF
        # lambdas = range(5433, 5436, step=5434.5/1e6)
        # prof = ones(length(lambdas))
        # lwavgrid = zeros(100)
        # rwavgrid = zeros(100)
        # allwavs = zeros(200)
        # allints = zeros(200)
        # GRASS.extrapolate_bisector(view(wav,:,1:1), view(bis,:,1:1))
        # GRASS.extrapolate_width(view(dep,:,1:1), view(wid,:,1:1))
        # plt.plot(wav[:,1], bis[:,1], c="tab:orange", ls="--")
        # GRASS.line_profile_cpu!(5434.535, lambdas, prof, wav[:,1] .- mean(wav[:,1]), dep[:,1], wid[:,1], lwavgrid, rwavgrid, allwavs, allints)

        # plt.plot(lambdas, prof, c="k", ls="--")
        # plt.show()
        # DEBUGGING STUFF

        # write input data to disk
        write_input_data(line_name, line_df.air_wavelength[1], fparams, wav, bis, dep, wid)
    end
    return nothing
end

function main()
    for name in line_info.name
        if name == "FeI_6173"
            println(">>> Processing " * name * "...")
            preprocess_line(name)
        end
    end
    return nothing
end

main()
