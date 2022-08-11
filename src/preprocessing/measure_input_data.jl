using CSV
using DataFrames
using Statistics
using GRASS

import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings

# read in table summarizing line info and spectra directories
line_info = CSV.read(GRASS.soldir * "line_info.csv", DataFrame)

# parse out columns
name = line_info.name
species = line_info.species
mass = line_info.mass
airwav = line_info.air_wav
geff = line_info.g_eff
height = line_info.height
lower = line_info.lower_level
upper = line_info.upper_level
spectra_dir = line_info.spectra_dir

# set input data absolute dir
data_dir = "/storage/home/mlp95/ford_dir/michael/lars_spectra/"

# plot average disk center spectra
for dir in unique(spectra_dir)
    df = GRASS.sort_spectrum_data(dir=data_dir * dir)
    df_dc = subset(df, :axis => x -> x .== "c")
    wavs, spec = GRASS.read_spectrum(df_dc.fpath[1] * df_dc.fname[1])
    wavs = mean(wavs, dims=2)
    spec = mean(spec, dims=2); spec ./= maximum(spec)
    plt.plot(wavs, spec); plt.show()
end
