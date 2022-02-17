# TODO check this line
#__precompile__
module GRASS

# parallelization modules
using CUDA; CUDA.allowscalar(false)
using Distributed
using SharedArrays

# import external modules
using CSV
using HDF5
using FITSIO
using Random
using StatsBase
using DataFrames
using Statistics
using SharedArrays
using StaticArrays
using Interpolations
import Glob.glob
import Dates.DateTime
import Polynomials: fit as pfit
import Base.Iterators: take, flatten, ProductIterator

# abbreviations for commonly used types
import Base: AbstractArray as AA
import Base: AbstractFloat as AF

# plots
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()

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
include("inputIO.jl")
include("bisectors.jl")

# rossiter mclaughlin
include("rm_effect.jl")

# star simulation
include("trim.jl")
include("synthesize.jl")
include("disk_sim.jl")

# processing spectra
include("velocities.jl")

# preprocessing of data
include("preprocessing/voigt.jl")
include("preprocessing/spectraIO.jl")
include("preprocessing/preprocessing.jl")

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
    synthesize_spectra(spec, disk; seed_rng=false, verbose=true)

Synthesize spectra given parameters in `spec` and `disk` instances.

# Arguments
- `spec::SpecParams`: SpecParams instance
- `disk::DiskParams`: DiskParams instance
"""
function synthesize_spectra(spec::SpecParams, disk::DiskParams; verbose::Bool=false,
                            use_gpu::Bool=false, seed_rng::Bool=false)
    # parse out dimensions for memory allocation
    N = disk.N
    Nλ = length(spec.lambdas)

    # allocate memory needed by both cpu & gpu implementations
    outspec = ones(Nλ, disk.Nt)

    # get number of calls to disk_sim needed
    indata_inds = unique(spec.data_inds)
    ncalls = length(indata_inds)

    # call appropriate simulation function on cpu or gpu
    if (use_gpu & CUDA.functional())
        # run the simulation and return
        for i in 1:ncalls
            # get temporary specparams with lines for this run
            spec_temp = SpecParams(spec, indata_inds[i])

            # load in the appropriate input data
            if verbose
                println("\t>>> " * spec.indata.dirs[indata_inds[i]])
            end
            soldata = SolarData(dir=spec.indata.dirs[indata_inds[i]]; spec.kwargs...)

            # run the simulation and multiply outspec by this spectrum
            disk_sim_gpu(spec_temp, disk, soldata, outspec, seed_rng=seed_rng, verbose=verbose)
        end
        return spec.lambdas, outspec
    else
        # allocate memory for synthesis
        prof = ones(Nλ)
        outspec_temp = zeros(Nλ, disk.Nt)

        # run the simulation (outspec modified in place)
        for i in eachindex(indata_inds)
            outspec_temp .= 0.0

            # get temporary specparams with lines for this run
            spec_temp = SpecParams(spec, indata_inds[i])

            # load in the appropriate input data
            if verbose
                println("\t>>> " * spec.indata.dirs[indata_inds[i]])
            end
            soldata = SolarData(dir=spec.indata.dirs[indata_inds[i]]; spec.kwargs...)

            # run the simulation and multiply outspec by this spectrum
            disk_sim(spec_temp, disk, soldata, prof, outspec_temp, seed_rng=seed_rng, verbose=verbose)
            outspec .*= outspec_temp
        end
        return spec.lambdas, outspec
    end
end

# precompile this function
precompile(synthesize_spectra, (SpecParams, DiskParams, Float64, Bool, Bool, Bool))

# export some stuff
export SpecParams, DiskParams, LineProperties, synthesize_spectra, calc_ccf, calc_rvs_from_ccf, calc_rms, use_gpu

end # module
