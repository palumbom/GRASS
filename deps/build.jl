using Tar, CodecZlib

moddir = abspath(joinpath(@__DIR__, ".."))
datdir = abspath(joinpath(moddir, "data/"))
soldir = abspath(joinpath(datdir, "input/"))

@assert isdir(moddir)
@assert isdir(datdir)
if !isdir(soldir)
    mkdir(soldir)
elseif !isempty(soldir)
    rm(soldir, recursive=true, force=true)
    mkdir(soldir)
end

download_url = "https://g-50c10.ffdaa9.e229.data.globus.org/GRASS_Input/input.tar.gz"
download_loc = abspath(joinpath(datdir, "input.tar.gz"))

download(download_url, download_loc)

Tar.extract(download_loc, soldir)
