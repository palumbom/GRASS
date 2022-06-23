# imports
using Pkg; Pkg.activate(".")
using CUDA
using GRASS
using Statistics

# plotting imports
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# define some functions
include(GRASS.moddir * "figures1/fig_functions.jl")

# get input data
inp = GRASS.InputData()
line_depths = GRASS.get_depth.(inp.lineprops)
