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

download_url = "https://zenodo.org/record/8271417/files/input_data.tar.gz?download=1"
download_loc = abspath(joinpath(datdir, "input.tar.gz"))

download(download_url, download_loc)

Tar.extract(download_loc, soldir)
