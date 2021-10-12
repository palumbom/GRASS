module GRASS

# parallelization modules
using CUDA
using Distributed
using SharedArrays

# import external modules
using CSV
using FITSIO
using Random
using StatsBase
using DataFrames
using Statistics
using SharedArrays
using Interpolations
import Glob.glob
import Dates.DateTime

# abbreviations for commonly used types
const AA = AbstractArray
const AF = AbstractFloat

# figure out if there is a GPU
const use_gpu = CUDA.functional()
if use_gpu
    const arr_type = CuArray
else
    const arr_type = Array
end

# configure directories
include("config.jl")

# ancillary functions + constants
include("utils.jl")
include("constants.jl")
include("interpolate.jl")

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
