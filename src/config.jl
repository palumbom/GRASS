# set absolute path to solar data
const moddir = abspath(joinpath(@__DIR__, ".."))
const datdir = abspath(joinpath(moddir, "data/"))
const soldir = abspath(joinpath(datdir, "input/"))
@assert isdir(moddir)
@assert isdir(soldir)
