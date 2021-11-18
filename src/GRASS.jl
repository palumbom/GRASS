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

# gpu implementation
include("gpu/gpu_utils.jl")
include("gpu/gpu_data.jl")
include("gpu/gpu_trim.jl")
include("gpu/gpu_sim.jl")
include("gpu/gpu_synthesis.jl")

"""
    synthesize_spectra(spec, disk; seed_rng=false, verbose=true, top=NaN)

Synthesize spectra given parameters in `spec` and `disk` instances.

# Arguments
- `spec::SpecParams`: SpecParams instance
- `disk::DiskParams`: DiskParams instance
"""
function synthesize_spectra(spec::SpecParams, disk::DiskParams;
                            top::Float64=NaN, seed_rng::Bool=false,
                            verbose::Bool=true, use_gpu::Bool=false)
    # parse out dimensions for memory allocation
    N = disk.N
    Nλ = length(spec.lambdas)

    # allocate memory for output
    outspec = zeros(Nλ, disk.Nt)

    # call appropriate simulation function on cpu or gpu
    if use_gpu
        # make sure there is actually a GPU
        @assert CUDA.functional()

        # run the simulation and return
        disk_sim_gpu(spec, disk, outspec)
        return spec.lambdas, outspec
    else
        # allocate memory for synthesis
        prof = ones(Nλ)

        # run the simulation (outspec modified in place)
        disk_sim(spec, disk, prof, outspec, seed_rng=seed_rng, verbose=verbose, top=top)
        return spec.lambdas, outspec
    end
end

# export some stuff
export SpecParams, DiskParams, synthesize_spectra, calc_ccf, calc_rvs_from_ccf, calc_rms, use_gpu

end # module
