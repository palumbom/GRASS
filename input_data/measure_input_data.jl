using Pkg; Pkg.activate(".")
using CSV
using HDF5
using DataFrames
using Statistics
using GRASS

import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
# mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# set LARS spectra absolute dir
data_dir = "/storage/group/ebf11/default/mlp95/lars_spectra/"

# read in table summarizing line info and spectra directories
line_info = CSV.read(GRASS.soldir * "line_info.csv", DataFrame)

# parse out columns
# name = line_info.name
# species = line_info.species
# mass = line_info.mass
# airwav = line_info.air_wav
# geff = line_info.g_eff
# height = line_info.height
# lower = line_info.lower_level
# upper = line_info.upper_level
# spectra_dir = line_info.spectra_dir

function write_line_params(line_name::String)
    # find row with line info
    line_df = subset(line_info, :name => x -> x .== line_name)

    # read in the IAG data and isolate the line
    iag_wavs, iag_flux = read_iag(isolate=true, airwav=line_df.air_wav[1])
    iag_depth = 1.0 - minimum(iag_flux)

    # get the filename
    line_dir = GRASS.soldir * line_name * "/"
    prop_file =  line_dir * line_name * "_line_properties.h5"

    # write the line properties file
    println(">>> Writing " * line_name * "_line_properties.h5")
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

function preprocess_line(line_name::String)
    # create subdirectory structure if it doesn't already exist
    if !isdir(GRASS.soldir * line_name)
        mkdir(GRASS.soldir * line_name)
    end



    return
end

function main()

    return
end
