# get Pkg manager stuff handled
using Pkg

# install PyPlot in global project (Khaled: uncomment on first run)
# Pkg.add("PyPlot")
using PyPlot

# Khaled: switch commenting on the next two lines (and change the path)
# Pkg.activate("PATH/TO/GRASS")
Pkg.activate(".")

# Khaled: for your first run uncomment the next two lines (and change the path)
# Khaled: you can comment these back out on subsequent runs
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
outdir = "/mnt/ceph/users/mpalumbo/data_for_khaled/"

# get the idx to run
template_idx = tryparse(Int, ARGS[1])

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

    # set some stuff for the ccf
    mask_type = EchelleCCFs.GaussianCCFMask # can also do top hot
    mask_width = (GRASS.c_ms/resolution) 
    Δv_step = 100.0 
    Δv_max = 32e3

    spec = SpecParams(lines=lines, depths=depths, variability=variability,
                      templates=templates, resolution=resolution, 
                      oversampling=4.0)

    # actually synthesize the spectra
    use_gpu = true # flip this if you don't have a GPU (it will be slow!)
    wavs, flux = synthesize_spectra(spec, disk, verbose=true, use_gpu=use_gpu)

    # measure bisectors (this uses linear interp)
    # (maybe one day I'll work on a smarter way)
    bis, int = GRASS.calc_bisector(wavs, flux, top=0.99, nflux=100)

    # calculate ccf 
    v_grid, ccf = calc_ccf(wavs, flux, lines, depths, resolution, 
                           mask_width=mask_width, mask_type=mask_type,
                           Δv_step=Δv_step, Δv_max=Δv_max)
    
    # measure velocities from ccf
    rvs, sigs = calc_rvs_from_ccf(v_grid, ccf, frac_of_width_to_fit=0.5)

    # write noiseless, full res spectra to disk
    fname = joinpath(outdir, lname[i] * "_noiseless.h5")
    h5open(fname, "w") do file
        write(file, "wavs", wavs, "flux", flux, "bis", bis, "int", int, "rvs", rvs) 
    end

    # convolve down the resolution
    new_res = 1.2e5
    oversampling = 4.0 # number of pixels per LSF FWHM (should be float not int)
    wavs_new, flux_new = GRASS.convolve_gauss(wavs, flux, new_res=new_res, 
                                              oversampling=oversampling)

    # make a copy to add noise to 
    flux_noisy = deepcopy(flux_new)
    snr = 500.0 # should be a float not int!
    GRASS.add_noise!(flux_noisy, snr)

    # measure bisectors
    bis_noisy, int_noisy = GRASS.calc_bisector(wavs_new, flux_new, top=0.99, nflux=100)

    # calculate ccf 
    v_grid_noisy, ccf_noisy = calc_ccf(wavs_new, flux_noisy, lines, depths, resolution, 
                                       mask_width=mask_width, mask_type=mask_type,
                                       Δv_step=Δv_step, Δv_max=Δv_max)
    
    # measure velocities from ccf
    rvs_noisy, sigs_noisy = calc_rvs_from_ccf(v_grid_noisy, ccf_noisy, frac_of_width_to_fit=0.5)

    # write the noisy data to disk
    fname = joinpath(outdir, lname[i] * "_noisy.h5")
    h5open(fname, "w") do file
        write(file, "wavs", wavs_new, "flux", flux_noisy, "bis", bis_noisy, "int", int_noisy, "rvs", rvs_noisy) 
    end

    # make a plot
    plot = false  
    if plot 
        plt.plot(wavs, flux[:,1], c="tab:blue", label="Noiseless")
        plt.plot(bis[:,1], int[:,1], c="tab:blue", ls="--")

        plt.plot(wavs_new, flux_noisy[:,1], c="tab:orange", label="SNR = " * string(snr))
        plt.plot(bis_noisy[:,1], int_noisy[:,1], c="tab:orange", ls="--")

        plt.xlabel("Wavelength (Å)")
        plt.ylabel("Flux")
        plt.legend()
        plt.show()
    end
end
