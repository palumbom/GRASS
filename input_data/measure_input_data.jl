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

# output directories for misc stuff
const grassdir, plotdir, datadir = GRASS.check_plot_dirs()
if !isdir(plotdir * "spectra_fits")
    mkdir(plotdir * "spectra_fits")
end

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

function find_wing_index(val, arr; min=argmin(arr))
    lidx = min - findfirst(x -> x .>= val, reverse(arr[1:min]))
    ridx = findfirst(x -> x .>= val, arr[min:end]) + min
    return lidx, ridx
end


function fit_line_wings(wavs_iso, flux_iso)
    # get indices and values for minimum, depth, and bottom
    min = argmin(flux_iso)
    bot = flux_iso[min]
    depth = 1.0 - bot

    # get wing indices for various percentage depths into line
    lidx50, ridx50 = find_wing_index(0.5 * depth + bot, flux_iso, min=min)
    lidx60, ridx60 = find_wing_index(0.6 * depth + bot, flux_iso, min=min)
    lidx70, ridx70 = find_wing_index(0.7 * depth + bot, flux_iso, min=min)
    lidx80, ridx80 = find_wing_index(0.8 * depth + bot, flux_iso, min=min)
    lidx90, ridx90 = find_wing_index(0.9 * depth + bot, flux_iso, min=min)

    # isolate the line wings and mask area around line core for fitting
    Δbot = 2
    core = min-Δbot:min+Δbot
    lwing = lidx90:lidx50
    rwing = ridx50:ridx90
    wavs_fit = vcat(wavs_iso[lwing], wavs_iso[core], wavs_iso[rwing])
    flux_fit = vcat(flux_iso[lwing], flux_iso[core], flux_iso[rwing])

    # set boundary conditions and initial guess
    # GOOD FOR FeI 5434 + others
    lb = [0.0, wavs_iso[min], 0.0, 0.0]
    ub = [1.0, wavs_iso[min], 0.5, 0.5]
    p0 = [1.0 - depth, wavs_iso[min], 0.02, 0.01]
    # GOOD FOR FeI 5434 + others

    # perform the fit
    fit = curve_fit(GRASS.fit_voigt, wavs_fit, flux_fit, p0, lower=lb, upper=ub)
    @show fit.param
    return fit
end

function replace_line_wings(fit, wavst, fluxt, min, val; debug=false)
    # get line model for all wavelengths in original spectrum
    flux_new = GRASS.fit_voigt(wavst, fit.param)

    # do a quick "normalization"
    flux_new ./= maximum(flux_new)

    # find indices
    idxl, idxr = find_wing_index(val, fluxt, min=min)

    # replace wings with model
    fluxt[1:idxl] .= flux_new[1:idxl]
    fluxt[idxr:end] .= flux_new[idxr:end]
    return nothing
end

function preprocess_line(line_name::String; verbose::Bool=true, debug::Bool=false)
    # create subdirectory structure if it doesn't already exist
    if !isdir(GRASS.soldir * line_name)
        mkdir(GRASS.soldir * line_name)
    end

    # find row with line info and write the line_params file
    line_df = subset(line_info, :name => x -> x .== line_name)
    if !debug
        write_line_params(line_df)
    end

    # find all the spectra files associated with this line
    fits_files = Glob.glob("*.fits", data_dir * line_df.spectra_dir[1] * "/")

    # read in the spectrum and bin into 15-second bins
    for i in eachindex(fits_files)
        # debugging block
        if debug && i > 1
            break
        end

        # print the filename
        if verbose
            println("\t >>> " * splitdir(fits_files[i])[end])
        end

        # get spec parameters
        fparams = GRASS.extract_line_params(fits_files[i])
        wavs, flux = GRASS.bin_spectrum(GRASS.read_spectrum(fits_files[i])...)

        # normalize the spectra
        flux ./= maximum(flux, dims=1)

        # allocate memory for input data
        wav = zeros(100, size(wavs,2))
        bis = zeros(100, size(wavs,2))
        wid = zeros(100, size(wavs,2))
        dep = zeros(100, size(wavs,2))

        # loop over epochs in spectrum file
        for t in 1:size(wavs, 2)
            # debugging block
            if debug && t > 1
                break
            end

            # get view of this time slice
            wavst = view(wavs, :, t)
            fluxt = view(flux, :, t)
            if debug
                fig, (ax1, ax2) = plt.subplots(1,2, figsize=(9,6))
                ax1.plot(wavst, fluxt, c="k", label="raw spec")
            end

            # refine the location of the minimum
            idx = findfirst(x -> x .>= line_df.air_wavelength[1], wavst)
            min = argmin(fluxt[idx-50:idx+50]) + idx - 50
            bot = fluxt[min]
            depth = 1.0 - bot

            # isolate the line
            idx1, idx2 = find_wing_index(0.95, fluxt, min=min)
            wavs_iso = copy(view(wavst, idx1:idx2))
            flux_iso = copy(view(fluxt, idx1:idx2))

            # fit the line wings
            fit = fit_line_wings(wavs_iso, flux_iso)

            # replace the line wings above val% continuum
            val = 0.9 #* depth + bot
            if debug
                idxl, idxr = find_wing_index(val, fluxt, min=min)
                ax1.axhline(fluxt[idxl], c="k", ls="--", alpha=0.5)
                ax1.axhline(fluxt[idxr], c="k", ls="--", alpha=0.5)
                ax1.plot(wavst, GRASS.fit_voigt(wavst, fit.param), c="tab:purple", label="model")
            end
            replace_line_wings(fit, wavst, fluxt, min, val, debug=debug)

            # TODO REVIEW BISECTOR CODE
            # measure the bisector and width function
            wav[:,t], bis[:,t] = GRASS.measure_bisector_interpolate(wavst, fluxt, top=0.99)
            dep[:,t], wid[:,t] = GRASS.measure_width_interpolate(wavst, fluxt, top=0.99)

            # debug plottings
            if debug
                ax1.plot(wavst, fluxt, c="tab:orange", ls="--", label="cleaned")
                ax1.plot(wav[:,t], bis[:,t], c="tab:blue")
                ax1.set_xlabel("Wavelength")
                ax1.set_ylabel("Normalized Intensity")
                ax1.legend()

                ax2.plot(dep[:,t], wid[:,t])
                ax2.set_xlabel("Normalized Intensity")
                ax2.set_ylabel("Width across line")
                fig.suptitle("\${\\rm " * replace(line_df.name[1], "_" => "\\ ") * "}\$")
                fig.savefig(plotdir * "spectra_fits/" * line_df.name[1] * ".pdf")
                plt.clf(); plt.close()
            end
        end

        # write input data to disk
        if !debug
            write_input_data(line_name, line_df.air_wavelength[1], fparams, wav, bis, dep, wid)
        end
    end
    return nothing
end

function main()
    for name in line_info.name
        name == "CI_5380" && continue
        name == "FeI_5382" && continue
        name == "NaI_5896" && continue
        name == "FeII_6149" && continue
        println(">>> Processing " * name)
        preprocess_line(name, debug=true)
    end
    return nothing
end

main()
