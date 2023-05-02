using Pkg; Pkg.activate(".")
using CSV
using CUDA
using GRASS
using Printf
using Revise
using DataFrames
using Statistics
using EchelleCCFs
using BenchmarkTools

using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")

# set up paramaters for spectrum
N = 132
Nt = 50
lines = [5250.6, 5434.5, 6173.4]
depths = [0.8, 0.9, 0.6]
geffs = [0.0, 0.0, 0.0]
templates = ["FeI_5250.6", "FeI_5434", "FeI_6173"]
variability = repeat([true], length(lines))
resolution = 7e5
seed_rng = true
use_gpu = true

disk = DiskParams(N=N, Nt=Nt)
spec1 = SpecParams(lines=lines, depths=depths, variability=variability,
                   geffs=geffs, templates=templates, resolution=resolution)
lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=seed_rng, verbose=true, use_gpu=use_gpu)

df = DataFrame("wave" => lambdas1, "flux" => dropdims(mean(outspec1, dims=2), dims=2))
CSV.write("/storage/work/m/mlp95/GRASS/test_spec.csv", df)
