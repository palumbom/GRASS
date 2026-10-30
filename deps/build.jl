using Tar, CodecZlib

moddir = abspath(joinpath(@__DIR__, ".."))
datdir = abspath(joinpath(moddir, "data/"))
soldir = abspath(joinpath(datdir, "input/"))

@assert isdir(moddir)
@assert isdir(datdir)

# kernel download helpers and env_flag (no SPICE dependency); download_kernels() uses the
# `datdir` defined above, and furnishing is deferred to GRASS.Eclipse's __init__ at runtime
include(joinpath(@__DIR__, "..", "src", "get_kernels.jl"))

# every build downloads the input data and tables again; GRASS_KEEP_DATA=true keeps the
# copies already on disk instead (CI restores them from a cache)
keep_data = env_flag("GRASS_KEEP_DATA")

# input data, extracted into an empty soldir
if !(keep_data && isdir(soldir) && !isempty(readdir(soldir)))
    rm(soldir, recursive=true, force=true)
    mkdir(soldir)

    download_url = "https://zenodo.org/record/8271417/files/input_data.tar.gz?download=1"
    download_loc = abspath(joinpath(datdir, "input.tar.gz"))

    download(download_url, download_loc)

    Tar.extract(download_loc, soldir)
end

files = [
    "LD_coeff_300.csv",
    "LD_coeff_SSD.csv",
    "LD_coeff_HD.csv",
    "sunspots.csv"
]

base_url = "https://zenodo.org/records/19372024/files/"

for f in files
    url = base_url * f * "?download=1"
    out = joinpath(datdir, f)

    keep_data && isfile(out) && continue

    download(url, out)

    # # Optional: auto-extract if it's a tar.gz
    # if endswith(f, ".tar.gz")
    #     Tar.extract(out, datadir)
    # end
end

# download the SPICE kernels
download_kernels()
