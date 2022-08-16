# imports
using CSV
using DataFrames
using Statistics
using GRASS

# plotting imports
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# set LARS spectra absolute dir and read line info file
const data_dir = "/storage/group/ebf11/default/mlp95/lars_spectra/"

# plot average disk center spectra
function plot_input_spectra(line_name::String, line_info::DataFrame; highlight::Bool=true)
    # find row with line info
    line_df = subset(line_info, :name => x -> x .== line_name)

    # get subset'd dataframes
    df = GRASS.sort_spectrum_data(dir=data_dir * line_df.spectra_dir[1] * "/")
    df_dc = subset(df, :axis => x -> x .== "c")
    lines = subset(line_info, :spectra_dir => x -> x .== splitdir(line_df.spectra_dir[1])[end])

    # make sure we are sorted on wavelength
    sort!(lines, :air_wavelength)

    # pull out line names and wavelengths
    names = lines.name
    airwavs = lines.air_wavelength

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
        idx = findfirst(x-> x .>= airwavs[i], wavs2)
        min = argmin(spec2[idx-50:idx+50]) + idx - 50

        # ax1.axvline(wavs2[min] ymin=0)
        if highlight
            if replace(line_name, "_" => " ") == replace(names[i], "_" => " ")
                ax1.annotate(("\${\\rm " * replace(names[i], "_" => "\\ ") * "}\$"),
                             (wavs2[min], spec2[min] - 0.025), (wavs2[min], spec2[min] - 0.15),
                             c="red", arrowprops=arrowprops, horizontalalignment="center",
                             fontsize=12)
            else
                ax1.annotate(("\${\\rm " * replace(names[i], "_" => "\\ ") * "}\$"),
                             (wavs2[min], spec2[min] - 0.025), (wavs2[min], spec2[min] - 0.15),
                             arrowprops=arrowprops, horizontalalignment="center",
                             fontsize=12)
            end
        else
            ax1.annotate(("\${\\rm " * replace(names[i], "_" => "\\ ") * "}\$"),
                         (wavs2[min], spec2[min] - 0.025), (wavs2[min], spec2[min] - 0.15),
                         arrowprops=arrowprops, horizontalalignment="center",
                         fontsize=12)
        end
    end

    # set labels and save the fig
    ax1.set_ylim(-0.1, 1.1)
    ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
    ax1.set_ylabel(L"{\rm Normalized\ Intensity}")
    fig.savefig(plotdir * "input_plots/" * line_name * "_spectrum.pdf")
    plt.clf(); plt.close()
    return nothing
end

function plot_input_data(line_name::String, line_info::DataFrame)
    # find row with line info
    line_df = subset(line_info, :name => x -> x .== line_name)

    # get directory of input data
    indir = GRASS.soldir * line_name * "/"
    soldata = GRASS.SolarData(dir=indir)

    key = (:c, :mu10)
    bis = soldata.bis[key]
    wav = soldata.wav[key]
    dep = soldata.dep[key]
    wid = soldata.wid[key]

    # create colormap and colorbar
    ncurves = 25
    nstp = round(Int, size(bis,2) / ncurves)
    iter = 1:nstp:size(bis,2)
    iter = 1:2:30
    cmap = plt.get_cmap("Blues", length(iter))
    cols = [cmap(1*i/length(iter)) for i in 1:length(iter)]
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=15 * iter[end]/60.0))

    # plot the bisectors
    fig = plt.figure()
    ax1 = fig.add_subplot()
    for (i, t) in enumerate(iter)
        ax1.plot(wav[:,t][2:end].*1000, bis[:,t][2:end], c=cols[i], alpha=0.75)
    end

    # shade upper region
    xlims = ax1.get_xlim()
    ax1.fill_between(range(xlims..., step=0.1), 0.8, 1.0, hatch="/",
                     fc="black", ec="white", alpha=0.15, zorder=0)

    # set axis limits + labels
    ax1.set_xlim(xlims...)
    ax1.set_ylim(0.1, 1.0)
    ax1.set_xlabel(L"{\rm Relative\ Wavelength\ (m\AA)}")
    ax1.set_ylabel(L"{\rm Normalized\ Intensity}")

    # set color bar
    cb = fig.colorbar(sm)
    cb.set_label(L"{\rm Time\ from\ first\ observation\ (min)}")

    # save the plot
    fig.savefig(plotdir * "input_plots/" * line_name * "_bisector.pdf")
    plt.clf(); plt.close()

    # plot the widths
    fig = plt.figure()
    ax1 = fig.add_subplot()
    for (i, t) in enumerate(iter)
        ax1.plot(dep[:,t], wid[:,t], c=cols[i], alpha=0.75)
    end

    # set axis limits + labels
    ax1.set_xlabel(L"{\rm Normalized\ Intensity}")
    ax1.set_ylabel(L"{\rm Width\ across\ line\ (\AA)}")

    # set color bar
    cb = fig.colorbar(sm)
    cb.set_label(L"{\rm Time\ from\ first\ observation\ (min)}")

    # save the plot
    fig.savefig(plotdir * "input_plots/" * line_name * "_width.pdf")
    plt.clf(); plt.close()
    return
end

function main()
    # make sure directory tree is set up
    if !isdir(plotdir * "input_plots/")
        mkdir(plotdir * "input_plots/")
    end

    # get summary of present lines
    line_info = CSV.read(GRASS.soldir * "line_info.csv", DataFrame)

    # loop over lines in line list
    for name in line_info.name
        # plot the spectrum
        plot_input_spectra(name, line_info; highlight=true)

        # plot the input data if it exists
        indir = GRASS.soldir * name * "/"
        if isdir(indir) && !isempty(indir)
            plot_input_data(name, line_info)
        end
    end
end

if (run | plot)
    main()
end

