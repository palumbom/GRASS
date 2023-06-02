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
using Distributions
using Interpolations
using PrecompileTools

# import specific methods
import Glob.glob
import Dates.DateTime
import Polynomials: fit as pfit, coeffs
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
include("preprocessing/conv_blueshift.jl")

# simulating observations
include("observing/convolutions.jl")
include("observing/signaltonoise.jl")
include("observing/ObservationPlan.jl")

# gpu implementation
include("gpu/gpu_utils.jl")
include("gpu/gpu_allocs.jl")
include("gpu/gpu_physics.jl")
include("gpu/gpu_data.jl")
include("gpu/gpu_trim.jl")
include("gpu/gpu_sim.jl")
include("gpu/gpu_synthesis.jl")

# functions for plotting figures
include("fig_functions.jl")
include("iag_utils.jl")


function generate_tloop!(tloop::AA{Int,2}, grid::StepRangeLen, soldata::SolarData{T}) where T<:AF
    # make sure dimensions are correct
    @assert size(tloop) == (length(grid), length(grid))

    # get the value of mu and ax codes
    disc_ax = parse_ax_string.(getindex.(keys(soldata.len),1))
    disc_mu = parse_mu_string.(getindex.(keys(soldata.len),2))

    # get indices to sort by mus
    inds_mu = sortperm(disc_mu)
    disc_mu .= disc_mu[inds_mu]
    disc_ax .= disc_ax[inds_mu]

    # get indices to sort by axis within mu sort
    for mu_val in unique(disc_mu)
        inds1 = (disc_mu .== mu_val)
        inds2 = sortperm(disc_ax[inds1])

        disc_mu[inds1] .= disc_mu[inds1][inds2]
        disc_ax[inds1] .= disc_ax[inds1][inds2]
    end

    for i in eachindex(grid)
        for j in eachindex(grid)
            # get positiosns
            x = grid[i]
            y = grid[j]

            # move to next iteration if off grid
            (x^2 + y^2) > one(T) && continue

            # get input data for place on disk
            key = get_key_for_pos(x, y, disc_mu, disc_ax)
            len = soldata.len[key]

            # generate random index
            tloop[i,j] = floor(Int, rand() * len) + 1
        end
    end
    return nothing
end

"""
    synthesize_spectra(spec, disk; seed_rng=false, verbose=true, top=NaN)

Synthesize spectra given parameters in `spec` and `disk` instances.

# Arguments
- `spec::SpecParams`: SpecParams instance
- `disk::DiskParams`: DiskParams instance
"""
function synthesize_spectra(spec::SpecParams{T}, disk::DiskParams{T};
                            seed_rng::Bool=false, verbose::Bool=true,
                            use_gpu::Bool=false, precision::DataType=Float64,
                            skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # parse out dimensions for memory allocation
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # allocate memory needed by both cpu & gpu implementations
    tloop_init = zeros(Int, N, N)
    outspec = ones(Nλ, Nt)

    # get number of calls to disk_sim needed
    templates = unique(spec.templates)

    # get grid
    grid = make_grid(N=disk.N)

    # call appropriate simulation function on cpu or gpu
    if use_gpu
        # make sure there is actually a GPU to use
        @assert CUDA.functional()

        # warn user if precision is single
        if precision <: Float32
           @warn "Single-precision GPU implementation produces large flux and velocity errors!"
        end

        # pre-allocate memory for gpu
        gpu_allocs = GPUAllocs(spec, disk, grid, precision=precision)

        # run the simulation and return
        for (idx, file) in enumerate(templates)
            # re-seed the rng
            if seed_rng
                Random.seed!(42)
            end

            # get temporary specparams with lines for this run
            spec_temp = SpecParams(spec, file)

            # load in the appropriate input data
            if verbose
                println("\t>>> " * splitdir(file)[end])
            end
            soldata = SolarData(fname=file)

            # generate or copy tloop
            if (idx > 1) && in_same_group(templates[idx - 1], templates[idx])
                @cusync CUDA.copyto!(gpu_allocs.tloop, tloop_init)
            else
                generate_tloop!(tloop_init, grid, soldata)
                @cusync CUDA.copyto!(gpu_allocs.tloop, tloop_init)
            end

            # run the simulation and multiply outspec by this spectrum
            disk_sim_gpu(spec_temp, disk, soldata, gpu_allocs, outspec, verbose=verbose,
                         skip_times=skip_times, precision=precision)
        end
        return spec.lambdas, outspec
    else
        # allocate memory for CPU
        prof = ones(Nλ)
        outspec_temp = zeros(Nλ, Nt)
        tloop = zeros(Int, N, N)

        # run the simulation (outspec modified in place)
        for (idx, file) in enumerate(templates)
            # re-seed the rng
            if seed_rng
                Random.seed!(42)
            end

            # get temporary specparams with lines for this run
            spec_temp = SpecParams(spec, file)

            # load in the appropriate input data
            if verbose
                println("\t>>> " * splitdir(file)[end])
            end
            soldata = SolarData(fname=file)

            # generate or copy tloop
            if (idx > 1) && in_same_group(templates[idx - 1], templates[idx])
                tloop .= tloop_init
            else
                generate_tloop!(tloop_init, grid, soldata)
                tloop .= tloop_init
            end

            # re-set array to 0s
            outspec_temp .= 0.0

            # run the simulation and multiply outspec by this spectrum
            disk_sim(spec_temp, disk, soldata, prof, outspec_temp, tloop,
                     skip_times=skip_times, verbose=verbose)
            outspec .*= outspec_temp
        end
        return spec.lambdas, outspec
    end
end

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
