using Tar
using Glob
using CodecZlib

moddir = abspath(joinpath(@__DIR__, ".."))
datadir = abspath(joinpath(moddir, "data"))
@assert isdir(datadir)

Tar.create(joinpath(datadir, "input"), joinpath(datadir, "input.tar.gz"))
