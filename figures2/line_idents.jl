using CSV
using DataFrames
using Statistics
using GRASS

import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

grassdir, plotdir, datadir = check_plot_dirs()

# read in table summarizing line info and spectra directories
line_info = CSV.read(GRASS.soldir * "line_info.csv", DataFrame)

# parse out columns
name = line_info.name
species = line_info.species
mass = line_info.mass
airwav = line_info.air_wav
geff = line_info.g_eff
height = line_info.height
lower = line_info.lower_level
upper = line_info.upper_level
spectra_dir = line_info.spectra_dir

# set input data absolute dir
data_dir = "/storage/home/mlp95/ford_dir/michael/lars_spectra/"

# plot average disk center spectra
function plot_line_idents(dir)
    # get subset'd dataframes
    df = GRASS.sort_spectrum_data(dir=dir)
    df_dc = subset(df, :axis => x -> x .== "c")
    lines = subset(line_info, :spectra_dir => x -> x .== splitdir(dir)[end])

    # make sure we are sorted on wavelength
    sort!(lines, :air_wav)

    # pull out line names and wavelengths
    names = lines.species
    waves = lines.air_wav

    # read in the spectra
    wavs, spec = GRASS.read_spectrum(df_dc.fpath[1] * df_dc.fname[1])
    wavs = dropdims(mean(wavs, dims=2), dims=2)
    spec = dropdims(mean(spec, dims=2), dims=2); spec ./= maximum(spec)

    # smooth the spectrum
    wavs = GRASS.moving_average(wavs, 5)
    spec = GRASS.moving_average(spec, 5)

    # interpolate up to high res for nice plotting of minima
    itp = GRASS.linear_interp(wavs, spec)
    wavs2 = range(minimum(wavs), maximum(wavs), step=minimum(diff(wavs)/5))
    spec2 = itp.(wavs2)

    # plot the spectra
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.plot(wavs2, spec2, "k", lw=1.5)

    # set arrow props
    arrowprops = Dict("facecolor"=>"black", "lw"=>1.5, "arrowstyle"=>"-")

    # annotate the lines
    for i in eachindex(names)
        # find the approx depth of the line
        idx = findfirst(x-> x .>= waves[i], wavs2)
        min = argmin(spec2[idx-50:idx+50]) + idx - 50

        # ax1.axvline(wavs2[min] ymin=0)
        ax1.annotate(names[i], (wavs2[min], spec2[min] - 0.025), (wavs2[min], spec2[min] - 0.15),
                     arrowprops=arrowprops, horizontalalignment="center")
    end

    # set labels and save the fig
    ax1.set_ylim(-0.1, 1.1)
    ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
    ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
    fig.savefig(plotdir * strip(splitdir(dir)[end], '/') * ".pdf")
    plt.clf(); plt.close()
    return nothing
end

plot_line_idents.(unique(data_dir .* spectra_dir));
