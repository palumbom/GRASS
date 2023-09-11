using Pkg; Pkg.activate(".")
using CUDA
using GRASS
using PyCall
using Printf
using Revise
using Statistics
using EchelleCCFs
using BenchmarkTools

# # plotting
# using LaTeXStrings
# import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
# mpl.style.use(GRASS.moddir * "fig.mplstyle")

# get astropy convolution
astroconv = pyimport("astropy.convolution")

# get synthetic spectrum
N = 132
Nt = 5
lines = [5434.5232]
depths = [0.9]
geffs = [0.0]
templates = ["FeI_5576"]
variability = repeat([true], length(lines))
resolution = 7e5
seed_rng = true

disk = DiskParams(N=N, Nt=Nt)
spec1 = SpecParams(lines=lines, depths=depths, variability=variability,
                   geffs=geffs, templates=templates, resolution=resolution)
lambdas1, outspec1 = synthesize_spectra(spec1, disk, seed_rng=false, verbose=true, use_gpu=true)
outspec1 = dropdims(mean(outspec1, dims=2), dims=2)

# get delta function
outspec2 = ones(length(outspec1))
outspec2[Int(length(lambdas1)/2)] = 0.1

# get the std in wavelength
std = 5434.5232 / 5e4 / 2.354

# get the std in pixels
pix_width = 5434.5232 / 7e5

kernel = astroconv.Gaussian1DKernel(stddev=std/pix_width)
convoluted1 = astroconv.convolve(outspec1, kernel, normalize_kernel=true, boundary="extend")
convoluted2 = astroconv.convolve(outspec2, kernel, normalize_kernel=true, boundary="extend")

# do my convolution
lambdas3, convoluted3 = GRASS.convolve_gauss(lambdas1, outspec1, new_res=5e4, oversampling=7e5/5e4)
lambdas4, convoluted4 = GRASS.convolve_gauss(lambdas1, outspec2, new_res=5e4, oversampling=7e5/5e4)

plt.plot(lambdas1, outspec1, label="og")
plt.plot(lambdas1, convoluted1, label="astropy")
plt.plot(lambdas3, convoluted3, label="grass")
plt.legend()
plt.show()
