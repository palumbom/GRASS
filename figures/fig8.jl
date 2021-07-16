using Pkg; Pkg.activate(".")
using Statistics
using GRASS
using EchelleCCFs
using FFTW

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

# function to get spectrum from input data
function spectrum_for_input(; mu::Symbol=:mu10, ax::Symbol=:c)
    # get all the input data
    soldata = GRASS.SolarData(extrapolate=true, contiguous_only=false)
    @assert haskey(soldata.bis, (ax, mu))

    # pull out necessary input and do line synthesis
    wav = soldata.wav[(ax, mu)]
    bis = soldata.bis[(ax, mu)]
    dep = soldata.dep[(ax, mu)]
    wid = soldata.wid[(ax, mu)]

    # set up arrays for synthesis
    lambdas = range(5434.5232 - 0.75, 5434.5232 + 0.75, step=5434.5232/7e5)
    lwavgrid = zeros(100)
    rwavgrid = zeros(100)
    allwavs  = zeros(200)
    allints  = zeros(200)
    flux = ones(length(lambdas), size(wav,2))
    prof = ones(length(lambdas))

    # synthesize the spectra
    for i in 1:size(wav,2)
        prof .= 1.0
        GRASS.line_from_bis!(5434.5232, lambdas, prof,
                             wav[:,i], dep[:,i], wid[:,i],
                             lwavgrid, rwavgrid, allwavs, allints)
        flux[:,i] .*= prof
    end
    return lambdas, flux
end

# function to get power spectrum for input data
function power_spec_for_input(; mu::Symbol=:mu10, ax::Symbol=:c)
    # first get the sepctrum
    lambdas, flux = spectrum_for_input(mu=mu, ax=ax)

    # now get the velocities
    v_grid, ccfs = calc_ccf(lambdas, flux, [5434.5232], [1-minimum(flux)], 7e5, normalize=true)
    rvs, sigs = calc_rvs_from_ccf(v_grid, ccfs)

    # do a fourier transform
    fourier = FFTW.fft(rvs) |> FFTW.fftshift
    freqs = FFTW.fftfreq(length(rvs), 1.0/15.0) |> FFTW.fftshift
    return freqs, abs.(fourier)
end

# now get the power spec for each disk position
mu = :mu10
ax = :c
freqs, power = power_spec_for_input(mu=mu, ax=ax)
plt.loglog(freqs, power, label="Input data")

# set up stuff for lines
N = 256
Nt = 1000
lines = [5434.5]
depths = [0.8]
variability = [true]
resolution = 700000.0
disk = DiskParams(N=N, Nt=Nt)
spec = SpecParams(lines=lines, depths=depths, variability=variability,
                  resolution=resolution, fixed_width=false,
                  fixed_bisector=false, extrapolate=true,
                  contiguous_only=false)

# synthesize spectra
println(">>> Synthesizing spectra...")
lambdas1, outspec1 = synthesize_spectra(spec, disk, seed_rng=false, verbose=true, top=NaN)

# calculate the ccf
println(">>> Calculating velocities...")
v_grid, ccf1 = calc_ccf(lambdas1, outspec1, spec, normalize=true)
rvs, sigs = calc_rvs_from_ccf(v_grid, ccf1)

# get frequencies to sample and then power
println(">>> Getting periodogram...")
fourier = FFTW.fft(rvs) |> FFTW.fftshift
freqs = FFTW.fftfreq(length(rvs), 1.0/15.0) |> FFTW.fftshift

# plot it
plt.loglog(freqs, abs.(fourier), label="Synthetic")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power (m/s)/Hz")
plt.legend()
plt.savefig("/Users/michael/Desktop/compare_ft.pdf")
plt.clf(); plt.close()
