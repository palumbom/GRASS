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

function __init__()
    # figure out if there is a gpu
    use_gpu = CUDA.functional()

    # initialize stuff for computations on GPU or CPU
    if use_gpu
        # define GPU function
        println(">>> Using GPU: " * CUDA.name(CUDA.device()))
        include(moddir * "src/gpu_functions.jl")

        # set array type to CuArray
        # const ArrayType = CuArray
        # time_loop = time_loop_gpu
        # line_loop = line_loop_gpu
        ArrayType = Array
        time_loop = time_loop_cpu
        line_loop = line_loop_cpu
        synth_func = line_profile_gpu!
    else
        # set array type to plain old array
        ArrayType = Array
        time_loop = time_loop_cpu
        line_loop = line_loop_cpu
        synth_func = line_profile_cpu!
    end
    return use_gpu, ArrayType, time_loop, line_loop, synth_func
end

use_gpu1, ArrayType1, time_loop1, line_loop1, synth_func1 = __init__()

const use_gpu = use_gpu1
const ArrayType = ArrayType1
const time_loop = time_loop1
const line_loop = line_loop1
const synth_func = synth_func1

# export some stuff
export SpecParams, DiskParams, synthesize_spectra, calc_ccf, calc_rvs_from_ccf, calc_rms

end # module
