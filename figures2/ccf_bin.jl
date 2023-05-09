using Pkg; Pkg.activate(".")
using CUDA
using JLD2
using GRASS
using Printf
using Revise
using FileIO
using DataFrames
using Statistics
using EchelleCCFs
using Polynomials
using Distributions
using BenchmarkTools
using HypothesisTests

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")

# make sure things are precompiled
function precompile()
    disk = DiskParams(N=132, Nt=5)
    spec1 = SpecParams(lines=[5250.6], depths=[0.5], templates=["FeI_5250.6"], resolution=7e5)
    lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=true, verbose=true, use_gpu=true)
    return nothing
end
precompile()

# get input data properties
lp = GRASS.LineProperties()
line_species = GRASS.get_species(lp)
rest_wavelengths = GRASS.get_rest_wavelength(lp)
line_depths = GRASS.get_depth(lp)
line_names = GRASS.get_name(lp)
line_titles = replace.(line_names, "_" => " ")
line_files = GRASS.get_file(lp)

# draw random line depths and centers
nlines = 10
lines = Float64[]
depths = Float64[]
templates = String[]
for i in eachindex(rest_wavelengths)
    # ltemp = rand(Uniform(minimum(rest_wavelengths), maximum(rest_wavelengths)), nlines)
    ltemp = rand(Uniform(5200, 5400), nlines)
    dtemp = rand(Normal(line_depths[i], 0.05), nlines)
    push!(lines, ltemp...)
    push!(depths, dtemp...)
    push!(templates, repeat([line_files[i]], nlines)...)
end

# synthesize a spectrum
N = 132
Nt = 1000
variability = trues(length(lines))
resolution = 7e5
seed_rng = true

disk = DiskParams(N=N, Nt=Nt)
spec1 = SpecParams(lines=lines, depths=depths, variability=variability, templates=templates, resolution=resolution)
lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=true, verbose=true, use_gpu=true)

# save it to a JLD
jldsave("spectra_for_bin.jld2", wavs=lambdas1, flux=outspec1)
