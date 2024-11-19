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

# configure directories
include("config.jl")

# get kernels for SPICE stuff
include("get_kernels.jl")

#set required body paramters as global variables 
#E,S,M radii (units:km)
earth_radius = bodvrd("EARTH", "RADII")[1]	
earth_radius_pole = bodvrd("EARTH", "RADII")[3]	
sun_radius = bodvrd("SUN","RADII")[1]
moon_radius = bodvrd("MOON", "RADII")[1] 

#collect LD info as global variables - (units: nm)
quad_ld_coeff_SSD = CSV.read(joinpath(datdir, "quad_ld_coeff_SSD.csv"), DataFrame)
quad_ld_coeff_300 = CSV.read(joinpath(datdir, "quad_ld_coeff_300.csv"), DataFrame)
quad_ld_coeff_HD = CSV.read(joinpath(datdir, "quad_ld_coeff_HD.csv"), DataFrame)
quad_ld_coeff_NL94 = CSV.read(joinpath(datdir, "quad_ld_coeff_NL94.csv"), DataFrame)

NL94_coeff = CSV.read(joinpath(datdir, "NL94_coeff.csv"), DataFrame)
lambda_nm = NL94_coeff.lambda_nm
a0 = NL94_coeff.a0
a1 = NL94_coeff.a1
a2 = NL94_coeff.a2
a3 = NL94_coeff.a3
a4 = NL94_coeff.a4
a5 = NL94_coeff.a5

extinction_coeff = DataFrame(CSV.File(joinpath(datdir, "NEID_three_extinction.csv")))

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
