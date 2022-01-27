# set absolute path to solar data
const moddir = abspath(joinpath(@__DIR__, ".."))
const soldir = abspath(joinpath(moddir, "input_data/"))
@assert isdir(moddir)
@assert isdir(soldir)
