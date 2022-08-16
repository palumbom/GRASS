using Pkg; Pkg.activate(".")
using CSV
using HDF5
using Glob
using DataFrames
using Statistics
using GRASS

import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# set LARS spectra absolute dir and read line info file
const data_dir = "/storage/group/ebf11/default/mlp95/lars_spectra/"
const line_info = CSV.read(GRASS.soldir * "line_info.csv", DataFrame)

# function to write line parameters file
function write_line_params(line_df::DataFrame)
    # get the filename
    line_dir = GRASS.soldir * line_df.name[1] * "/"
    prop_file =  line_dir * line_df.name[1] * "_line_properties.h5"

    # don't bother if the file already exists
    if isfile(prop_file)
        println("\t >>> " * splitdir(prop_file)[end] * " already exists...")
        return nothing
    end

    # read in the IAG data and isolate the line
    iag_wavs, iag_flux = read_iag(isolate=true, airwav=line_df.air_wav[1])
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

function write_input_data(line_name, air_wav, fparams, wav, bis, dep, wid)
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
        attr["wavelength"] = air_wav

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
    for f in fits_files
        # print the filename
        if verbose println("\t >>> Processing " * splitdir(f)[end]) end

        # DEBUGGING STUFF
        """
        if splitdir(f)[end] != "lars_l12_20180427-114647_clv5576_mu09_s.ns.chvtt.fits"
            continue
        end
        """

        # get spec parameters
        fparams = GRASS.extract_line_params(f)
        wavs, flux = GRASS.bin_spectrum(GRASS.read_spectrum(f)...)

        # fparams = GRASS.extract_line_params(fits_files[1])
        # wavs, flux = GRASS.bin_spectrum(GRASS.read_spectrum(fits_files[1])...)

        # normalize the spectra
        flux ./= maximum(flux, dims=1)

        # loop over epochs in data and measure input data
        wav = zeros(100, size(wavs,2))
        bis = zeros(100, size(wavs,2))
        wid = zeros(100, size(wavs,2))
        dep = zeros(100, size(wavs,2))
        for t in 1:size(wavs, 2)
            # refine the location of the minimum
            idx = findfirst(x -> x .>= line_df.air_wav[1], wavs[:,t])
            min = argmin(flux[idx-50:idx+50, t]) + idx - 50

            # isolate the line
            idx1 = findfirst(x -> x .>= wavs[min - 75, t], wavs[:, t])
            idx2 = findfirst(x -> x .>= wavs[min + 75, t], wavs[idx1:end, t]) + idx1
            wavs_iso = view(wavs, idx1:idx2, t)
            flux_iso = view(flux, idx1:idx2, t)

            # measure a bisector
            wav[:,t], bis[:,t] = GRASS.measure_bisector(wavs_iso, flux_iso, interpolate=false)
            dep[:,t], wid[:,t] = GRASS.measure_width_loop(wavs_iso, flux_iso)

            # TODO REVIEW BISECTOR CODE
            # wav[:,t], bis[:,t] = GRASS.calc_bisector(wavs_iso, flux_iso)
            # dep[:,t], wid[:,t] = GRASS.calc_width_function(wavs_iso, flux_iso)
        end

        # write input data to disk
        write_input_data(line_name, line_df.air_wav[1], fparams, wav, bis, dep, wid)
    end
    return nothing
end

function main()
    for name in line_info.name
        if name == "NiI_5578"
            println(">>> Processing " * name * "...")
            preprocess_line(name)
        end
    end
    return nothing
end

main()
