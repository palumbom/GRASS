using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using FFTW
@everywhere using GRASS
@everywhere using Statistics
@everywhere using EchelleCCFs
@everywhere using SharedArrays
using JLD2
using FileIO
using LaTeXStrings

# define rms loop function
@everywhere include(GRASS.moddir * "figures/fig_functions.jl")

# some global stuff
const N = 256
const Nt = 500
const Nloop = 200

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# function to get spectrum from input data
@everywhere function spectrum_for_input(; mu::Symbol=:mu10, ax::Symbol=:c)
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
@everywhere function power_spec_for_input(; mu::Symbol=:mu10, ax::Symbol=:c)
    # first get the sepctrum
    lambdas, flux = spectrum_for_input(mu=mu, ax=ax)

    # now get the velocities
    v_grid, ccfs = calc_ccf(lambdas, flux, [5434.5232], [1-minimum(flux)], 7e5, normalize=true)
    rvs, sigs = calc_rvs_from_ccf(v_grid, ccfs)

    # get the power spectrum and return
    return power_spectrum(15.0, rvs)
end

@everywhere function power_spectrum(period, signal)
    # do fourier transform and get frequencies
    fourier = FFTW.fft(signal)
    freqs = range(0.0, 1.0/(2.0 * period), length=length(signal) ÷ 2)

    # get power
    power = abs.(fourier).^2 ./ (freqs[2] - freqs[1])
    power = (2.0 / length(signal)) * power[1:length(signal) ÷ 2]
    # power = power[1:length(signal) ÷ 2]
    return freqs, power
end

function main()
    """
    # now get the power spec for each disk position
    mu = :mu10
    ax = :c
    freqs_dat, power_dat = power_spec_for_input(mu=mu, ax=ax)
    plt.loglog(freqs_dat, power_dat, label="Input data")
    """

    # set up spectrum parameters
    lines = [5434.5]
    depths = [0.8]
    variability = [true]
    resolution = 700000.0
    disk = DiskParams(N=N, Nt=Nt)
    spec = SpecParams(lines=lines, depths=depths, variability=variability,
                      resolution=resolution, fixed_width=false,
                      fixed_bisector=false, extrapolate=true,
                      contiguous_only=false)

    # allocate memory
    freqs = SharedArray{Float64}(Nt ÷ 2, Nloop)
    powers = SharedArray{Float64}(Nt ÷ 2, Nloop)

    # calculate power spectrum many times to build up stats
    @sync @distributed for i in 1:Nloop
        # synthesize spectra and compute power spectrum
        lambdas1, outspec1 = synthesize_spectra(spec, disk, seed_rng=false)
        v_grid, ccf1 = calc_ccf(lambdas1, outspec1, spec, normalize=true)
        rvs, sigs = calc_rvs_from_ccf(v_grid, ccf1)
        freqs_sim, power_sim = power_spectrum(15.0, rvs)

        # assign to shared array
        freqs[:, i] = freqs_sim
        powers[:, i] = power_sim
    end

    # convert to plain arrays
    freqs = convert(Array{Float64}, freqs)
    powers = convert(Array{Float64}, powers)

    # save the output
    outfile = datadir * "power_loop_" * string(Nloop) * ".jld2"
    save(outfile,
         "freqs", freqs,
         "powers", powers)
    return nothing
end

# run the simulation
if run
    main()
end

# plotting code block
if plot
    # plotting imports
    import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
    using PyCall; animation = pyimport("matplotlib.animation")
    mpl.style.use(GRASS.moddir * "figures/fig.mplstyle")

    # read in the data
    file = datadir * "power_loop_" * string(Nloop) * ".jld2"
    d = load(file)
    freqs = d["freqs"]
    powers = d["powers"]

    # compute mean and std
    avg_power = mean(2/500 .* powers, dims=2)
    std_power = std(2/500 .* powers, dims=2)

    # plot it
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.loglog(freqs[:,1], avg_power)
    # ax1.fill_between(freqs[:,1], log10.(avg_power .- std_power), log10.(avg_power .+ std_power), color="tab:blue", alpha=0.3)

    # set labels, etc.
    ax1.set_ylim(1.0, 1e6)
    ax1.set_xlabel(L"{\rm Frequency\ (Hz)}")
    ax1.set_ylabel(L"{\rm Power\ (arbitrary\ units)}")

    # save the figure
    fig.savefig(plotdir * "fig8.pdf")
    plt.clf(); plt.close()
    println(">>> Figure written to: " * plotdir * "fig8.pdf")
end
