using GRASS
using CUDA
using Test

# run the CPU tests
include("test_geometry.jl")
include("test_physics.jl")
include("test_ccf.jl")
include("test_interpolations.jl")   # TODO: obsoleted?
include("test_input.jl")
include("test_synthesis.jl")

# run the GPU tests if there is a GPU
if CUDA.functional()
    include("test_gpu.jl")
end
