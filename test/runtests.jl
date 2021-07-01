using GRASS
using Test

# run the tests
include("test_geometry.jl")
include("test_physics.jl")
include("test_ccf.jl")
include("test_interpolations.jl")
include("test_input.jl")  # TODO figure out this issue with automated testing
include("test_synthesis.jl")
