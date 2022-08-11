using Pkg; Pkg.activate(".")
using CSV
using DataFrames
using Statistics
using GRASS

import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
# mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# plotting functions
include(GRASS.moddir * "figures1/fig_functions.jl")
grassdir, plotdir, datadir = check_plot_dirs()

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
