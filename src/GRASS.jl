module GRASS

# parallelization modules
using CUDA; CUDA.allowscalar(false)
using Distributed
using SharedArrays

# import external modules
using SPICE
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
import NaNMath
import NaNMath: sum as nansum
import Polynomials: fit as pfit, coeffs

# abbreviations for commonly used types
import Base: AbstractArray as AA
import Base: AbstractFloat as AF

#set required body paramters as global variables 
#E,S,M radii (units:km)
include("get_kernels.jl")
earth_radius = bodvrd("EARTH", "RADII")[1]	
earth_radius_pole = bodvrd("EARTH", "RADII")[3]	
sun_radius = bodvrd("SUN","RADII")[1]
moon_radius = bodvrd("MOON", "RADII")[1] 

lambda_nm = [303.327, 310.843, 320.468, 329.897, 349.947, 365.875, 374.086, 390.915, 401.970, 416.319, 427.930, 443.885, 445.125, 457.345, 477.427, 492.905, 519.930, 541.76, 559.95, 579.88, 610.975, 640.97, 669.4, 700.875, 748.71, 811.76, 869.6, 948.85, 1046.6, 1098.95]
a0 = [0.08011, 0.08160, 0.08833, 0.09188, 0.11012, 0.12828, 0.12579, 0.12995, 0.12323, 0.12814, 0.14249, 0.16220, 0.15248, 0.16604, 0.19571, 0.20924, 0.23695, 0.26073, 0.26892, 0.28392, 0.30854, 0.33644, 0.34685, 0.37885, 0.40627, 0.42977, 0.45396, 0.47855, 0.49870, 0.51149]
a1 = [0.70695, 0.71609, 0.77285, 0.92459, 1.02168, 1.04969, 0.85402, 0.91836, 1.08648, 1.19947, 1.28796, 1.24893, 1.38517, 1.38544, 1.30551, 1.30798, 1.29927, 1.27428, 1.34319, 1.36896, 1.3662, 1.30590, 1.37539, 1.25553, 1.22842, 1.25182, 1.25101, 1.19813, 1.21429, 1.19354]
a2 = [0.4991, 0.69685, 0.65382, 0.19604, -0.10924, 0.17482, 0.54601, -0.07566, -0.43974, -0.84407, -1.19564, -0.92165, -1.49615, -1.52275, -1.25845, -1.20411, -1.28034, -1.30352, -1.58427, -1.75998, -1.83572, -1.79238, -2.04425, -1.70908, -1.67877, -1.85164, -2.02958, -1.86296, -2.06976, -2.00174]
a3 = [-0.31080, -0.87703, -1.04647, -0.39546, -0.00055, -1.13371, -1.15048, 0.19149, 0.45912, 1.07201, 1.68603, 0.89978, 1.99886, 2.00232, 1.50626, 1.21505, 1.37760, 1.47085, 1.91271, 2.22154, 2.33221, 2.4504, 2.70493, 2.19647, 2.05535, 2.31949, 2.7541, 2.36939, 2.80703, 2.66936]
a4 = [-0.02177, 0.47008, 0.72921, 0.23599, -0.08688, 1.23882, 0.88928, -0.28712, -0.32759, -0.79537, -1.36658, -0.50148, -1.48155, -1.45969, -1.05472, -0.67196, -0.85054, -0.96618, -1.3135, -1.56074, -1.63082, -1.89979, -1.9429, -1.59554, -1.39972, -1.59101, -2.02287, -1.64367, -2.05247, -1.94981]
a5 = [0.04642, -0.0876, -0.19775, -0.05303, 0.06487, -0.4599, -0.26462, 0.12298, 0.0985, 0.23982, 0.44572, 0.1122, 0.44119, 0.42864, 0.30570, 0.14381, 0.21706, 0.26384, 0.37295, 0.4463, 0.45959, 0.59943, 0.55999, 0.47378, 0.38845, 0.44155, 0.59338, 0.46056, 0.60221, 0.57715]

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
include("star_geometry.jl")
include("star_physics.jl")

# geometry for orbiting bodies
include("state_vectors.jl")
include("kepler_equation.jl")

# data read-in + calculations
include("inputIO.jl")
include("bisectors.jl")

# star simulation
include("trim.jl")
include("synthesize.jl")
include("disk_sim.jl")
include("eclipse_comp.jl")
include("disk_sim_eclipse.jl")
include("disk_precomps.jl")

# processing spectra
include("velocities.jl")

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
include("gpu/gpu_precomps_eclipse.jl")
include("gpu/gpu_trim.jl")
include("gpu/gpu_sim.jl")
include("gpu/gpu_sim_eclipse.jl")
include("gpu/gpu_synthesis.jl")
include("gpu/gpu_state_vectors.jl")

# functions for plotting figures
include("fig_functions.jl")
include("iag_utils.jl")

# include convenience functions for synthtesis
include("convenience.jl")
include("convenience_rossiter.jl")
include("convenience_eclipse.jl")

# export some stuff
export SpecParams, DiskParams, DiskParamsEclipse, LineProperties, SolarData, Planet,
       synthesize_spectra, simulate_rossiter, calc_ccf,
       calc_rvs_from_ccf, calc_rms, parse_args,
       check_plot_dirs, read_iag

end # module
