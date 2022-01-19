# set absolute path to solar data
const moddir = abspath(joinpath(@__DIR__, ".."))
const soldir = abspath(joinpath(moddir, "old_input_data/FeI_5434/"))
@assert isdir(moddir)
@assert isdir(soldir)
