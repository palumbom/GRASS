module GRASS # parent module

# parallelization
using CUDA; CUDA.allowscalar(false)

# import external packages
using CSV
using HDF5
using JLD2
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
using ProgressMeter
using ImageFiltering
using Interpolations
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

# include convenience functions for synthtesis
include("convenience.jl")
include("resolved.jl")

# export some stuff
export SpecParams, DiskParams, LineProperties, SolarData, synthesize_spectra,
       calc_ccf, calc_rvs_from_ccf, calc_rms, parse_args, check_plot_dirs,
       read_iag, calc_bisector, calc_bisector_inverse_slope, calc_bisector_span,
       calc_bisector_bottom, calc_bisector_curvature, moving_average#,
       #DiskParamsEclipse, simulate_rossiter, Planet

# module GRASSe # eclipse submodule
module Eclipse # eclipse submodule

# inherit from parent module
using CSV
using SPICE
using GRASS
using CUDA
using DataFrames
using Statistics
datdir = GRASS.datdir

import Base: AbstractArray as AA
import Base: AbstractFloat as AF

# get kernels for SPICE stuff
include("get_kernels.jl")

#set required body paramters as global variables 
#E,S,M radii (units:km)
earth_radius = bodvrd("EARTH", "RADII")[1]	
earth_radius_pole = bodvrd("EARTH", "RADII")[3]	
sun_radius = bodvrd("SUN","RADII")[1]
moon_radius = bodvrd("MOON", "RADII")[1] 

#collect LD info as global variables - (units: nm)
quad_ld_coeff_SSD = CSV.read(joinpath(datdir, "LD_coeff_SSD.csv"), DataFrame)
quad_ld_coeff_300 = CSV.read(joinpath(datdir, "LD_coeff_300.csv"), DataFrame)
quad_ld_coeff_HD = CSV.read(joinpath(datdir, "LD_coeff_HD.csv"), DataFrame)

spots_info = DataFrame(CSV.File(joinpath(datdir, "sunspots.csv")))

# structures 
include("structures/DiskParamsEclipse.jl")
include("structures/SynthWorkspaceEclipse.jl")
include("structures/GPUAllocsEclipse.jl")
include("structures/GPUSolarData.jl")

# eclipse stuff
include("synthesize_eclipse.jl")
include("eclipse_comp.jl")
include("disk_sim_eclipse.jl")
include("convenience_eclipse.jl")

# gpu implementation
include("gpu/gpu_physics_eclipse.jl")
include("gpu/gpu_precomps_eclipse.jl")
include("gpu/gpu_sim_eclipse.jl")

export synthesize_spectra_eclipse

end # submodule

end # parent module
