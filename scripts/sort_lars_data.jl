using Pkg; Pkg.activate(".")
using CSV
using Glob
using DataFrames
using Statistics
using GRASS

# set input data directory
lars_dir = abspath("/mnt/ceph/users/mpalumbo/LARS_spectra")

# set output directory
out_dir = joinpath(lars_dir, "all_fits")
if !isdir(out_dir)
    mkdir(out_dir)
end

# loop over files in directory
for dd in readdir(lars_dir, join=true)
    # skip the output dir and the tar archives
    if (dd == out_dir) | (!isdir(dd))
        continue
    end

    # now copy the file over
    l1_file = readdir(joinpath(dd, "l1", "l1data"), join=true)[1]
    cp(l1_file, joinpath(out_dir, basename(l1_file)), force=true)
end