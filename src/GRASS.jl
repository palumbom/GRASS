module GRASS

# import external modules
using CSV
using LsqFit
using FITSIO
using Random
# using Debugger
using StatsBase
using DataFrames
using Statistics
using SharedArrays
# using LinearAlgebra
using Interpolations
import Glob.glob
import Dates.DateTime

# abbreviations for commonly used types
const AA = AbstractArray
const AF = AbstractFloat

# configure directories
include("config.jl")

# ancillary functions + constants
include("utils.jl")
include("constants.jl")

# user-defined composite types
include("structures.jl")

# star geometry + thermal/RT physics
include("stargeometry.jl")
include("starphysics.jl")

# data read-in + calculations
include("lineIO.jl")
include("bisectors.jl")

# star simulation
include("trim.jl")
include("synthesize.jl")
include("disk_sim.jl")

# processing spectra
include("velocities.jl")

# simulating observations
include("observing/convolutions.jl")
include("observing/signaltonoise.jl")
include("observing/ObservationPlan.jl")

# export some stuff
export SpecParams, DiskParams, synthesize_spectra, calc_ccf, calc_rvs_from_ccf, calc_rms

end # module
