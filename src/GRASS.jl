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
include("ccfs/rv_uncertainty.jl")

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
using LinearAlgebra
using Random
const datdir = GRASS.datdir

import Base: AbstractArray as AA
import Base: AbstractFloat as AF

# get kernels for SPICE stuff (defines download_kernels() and furnsh_kernels())
include("get_kernels.jl")
include("earth_frame.jl")

# body radii (km) and limb-darkening / sunspot tables — declared here with concrete
# types for type-stable access, but populated at runtime in __init__ (not at precompile)
global earth_radius::Float64
global earth_radius_pole::Float64
global sun_radius::Float64
global moon_radius::Float64
global quad_ld_coeff_SSD::DataFrame
global quad_ld_coeff_300::DataFrame
global quad_ld_coeff_HD::DataFrame
global spots_info::DataFrame

# furnish SPICE kernels and load data tables at runtime, NOT at module-load/precompile time
function __init__()
    # download_kernels() is a no-op once Pkg.build has fetched the kernels; it guards the
    # kernels-missing-but-input-present case (config.jl's self-heal only checks data/input/)
    download_kernels()

    # the high-precision Earth PCK expires about three months after NAIF generates it
    if earth_pck_stale()
        try
            download_earth_pck!()
        catch err
            @warn "Could not refresh earth_latest_high_prec.bpc from NAIF; ITRF93 may lack coverage for recent epochs" exception=err
        end
    end
    furnsh_kernels()

    # E, S, M radii (units: km) — requires the kernel pool to be furnished above
    global earth_radius = bodvrd("EARTH", "RADII")[1]
    global earth_radius_pole = bodvrd("EARTH", "RADII")[3]
    global sun_radius = bodvrd("SUN", "RADII")[1]
    global moon_radius = bodvrd("MOON", "RADII")[1]

    # limb-darkening coefficients (wavelength in nm) and the sunspot table: lon and lat in
    # degrees as body-fixed IAU_SUN angles (Carrington-like), diameter_km, diameter_arcsec,
    # and contrast as the fractional intensity deficit (I_phot - I_spot) / I_phot
    global quad_ld_coeff_SSD = CSV.read(joinpath(datdir, "LD_coeff_SSD.csv"), DataFrame)
    global quad_ld_coeff_300 = CSV.read(joinpath(datdir, "LD_coeff_300.csv"), DataFrame)
    global quad_ld_coeff_HD = CSV.read(joinpath(datdir, "LD_coeff_HD.csv"), DataFrame)
    global spots_info = DataFrame(CSV.File(joinpath(datdir, "sunspots.csv")))

    return nothing
end

include("utils.jl")

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
include("gpu/gpu_synthesis_eclipse.jl")
include("gpu/gpu_sim_eclipse.jl")

# star geometry + thermal/RT physics
include("star_geometry.jl")
include("star_physics.jl")

export synthesize_spectra_eclipse

end # submodule

end # parent module
