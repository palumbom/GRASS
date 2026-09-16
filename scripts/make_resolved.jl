# get Pkg manager stuff handled
using Pkg

# install PyPlot in global project
# Pkg.add("PyPlot")
using PyPlot

# Pkg.activate("PATH/TO/GRASS")
Pkg.activate(".")

# Pkg.develop(path="PATH/TO/GRASS")
# Pkg.instantiate()

# import packages
using CSV
using HDF5
using GRASS
using Revise
using DataFrames
using EchelleCCFs

# output directory
outdir = "/mnt/ceph/users/mpalumbo/data_for_eduardo/"
!isdir(outdir) && mkdir(outdir)

# get the idx to run
if !isempty(ARGS)
    show_progress = false
    template_idx = tryparse(Int, ARGS[1])
else
    show_progress = true
    template_idx = 9
end

# get data for the template lines in GRASS library
lp = GRASS.LineProperties()
wavelength = lp.λrest
depth = lp.depth
dfile = lp.file
lname = GRASS.get_name(lp)

# set up parameters for synthetic spectrum
Nt = 24000 # number of 15-second time steps in simulation
disk = DiskParams(Nt=Nt)

# loop over lines in library
let i = template_idx
    # parameters for synthetic spectrum
    lines = [wavelength[i]]
    depths = [depth[i]]
    templates = [dfile[i]]
    variability = trues(length(lines)) # control whether lines dance
    resolution = 7e5 # don't set resolution here, rather convolve it down later

    # spec object
    spec = SpecParams(lines=lines, depths=depths, variability=variability,
                      templates=templates, resolution=resolution, 
                      oversampling=4.0)

    # set mu bins
    μ_bins = range(0.025, 0.975, step=0.05)

    # synthesize
    wavs, flux = GRASS.synthesize_spectra_resolved(μ_bins, spec, disk, verbose=true, 
                                                   use_gpu=true, show_progress=show_progress)

    # write noiseless, full res spectra to disk
    fname = joinpath(outdir, lname[i] * ".h5")
    h5open(fname, "w") do file
        write(file, "wavs", wavs, "flux", flux, "mu_bin_center", collect(μ_bins))
    end
end
