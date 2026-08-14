using GRASS
using CUDA
using Test

# run the CPU tests
include("test_geometry.jl")
include("test_physics.jl")
include("test_ccf.jl")
include("test_rv_uncertainty.jl")
include("test_interpolations.jl")
include("test_input.jl")
include("test_synthesis.jl")

# run the GPU tests if there is a GPU; CI has none, so make the skip visible
if CUDA.functional()
    include("test_gpu.jl")
else
    @info "No functional GPU detected: CPU/GPU parity tests SKIPPED (test_gpu.jl)"
end
