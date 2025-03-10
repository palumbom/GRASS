module GRASS

# parallelization modules
using CUDA; CUDA.allowscalar(false)
using Distributed
using SharedArrays

# import external modules
using CSV
using HDF5
using Dates
using LsqFit
using FITSIO
using Random
using Dierckx
using StatsBase
using DataFrames
using Statistics
using Polynomials
using LinearAlgebra
using Distributions
using ImageFiltering
using Interpolations
using PrecompileTools
using OrderedCollections

# import specific methods
import Glob.glob
import Dates.DateTime
import NaNMath: sum as nansum
import Polynomials: fit as pfit, coeffs

# abbreviations for commonly used types
import Base: AbstractArray as AA
import Base: AbstractFloat as AF

# configure directories
include("config.jl")

# ancillary functions + constants
include("utils.jl")
include("gpu/gpu_utils.jl")
include("constants.jl")
include("interpolate.jl")

# composite types
include("structures.jl")

# star geometry + thermal/RT physics
include("stargeometry.jl")
include("starphysics.jl")

# data read-in + calculations
include("inputIO.jl")
include("bisectors.jl")

# star simulation
include("trim.jl")
include("synthesize.jl")
include("disk_sim.jl")
include("disk_precomps.jl")

# processing spectra
include("ccfs/ccf.jl")

# preprocessing of data
include("preprocessing/voigt.jl")
include("preprocessing/spectraIO.jl")
include("preprocessing/preprocessing.jl")
include("preprocessing/conv_blueshift.jl")

# simulating observations
include("observing/convolutions.jl")
include("observing/signaltonoise.jl")
include("observing/ObservationPlan.jl")

# gpu implementation
include("gpu/gpu_physics.jl")
include("gpu/gpu_data.jl")
include("gpu/gpu_precomps.jl")
include("gpu/gpu_trim.jl")
include("gpu/gpu_sim.jl")
include("gpu/gpu_synthesis.jl")

# functions for plotting figures
include("fig_functions.jl")
include("iag_utils.jl")

# lastly include convenciece synthesis
include("convenience.jl")

# precompile stuff
# if !isempty(soldir)
#     @compile_workload begin
#         # params for spectra
#         N = 132
#         Nt = 2
#         lines = [5434.5]
#         depths = [0.5]
#         geffs = [0.0]
#         templates = ["FeI_5434"]
#         variability = repeat([true], length(lines))
#         resolution = 7e5

#         # make composite type instances
#         disk = DiskParams(N=N, Nt=Nt)
#         spec = SpecParams(lines=lines, depths=depths, variability=variability,
#                            geffs=geffs, templates=templates, resolution=resolution)

#         if CUDA.functional()
#             # TODO figure out what's going wrong with GPU precompilation caching
#             # lambdas1, outspec = synthesize_spectra(spec, disk, seed_rng=true, verbose=false, use_gpu=true)
#             lambdas1, outspec = synthesize_spectra(spec, disk, seed_rng=true, verbose=false, use_gpu=false)
#         else
#             lambdas1, outspec = synthesize_spectra(spec, disk, seed_rng=true, verbose=false, use_gpu=false)
#         end
#     end
# end

# export some stuff
export SpecParams, DiskParams, LineProperties, SolarData, synthesize_spectra,
       calc_ccf, calc_rvs_from_ccf, calc_rms, parse_args, check_plot_dirs,
       read_iag

end # module
