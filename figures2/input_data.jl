# imports
using CSV
using DataFrames
using Statistics
using GRASS
using PyCall

# plotting imports
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
mpl.style.use(joinpath(GRASS.moddir, "fig.mplstyle"))
adjust_text = pyimport("adjustText")

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
    fname = GRASS.soldir * line_name * ".h5"
    soldata = GRASS.SolarData(fname)

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

# function main()
    # make sure directory tree is set up
    if !isdir(plotdir * "input_plots/")
        mkdir(plotdir * "input_plots/")
    end

    # # get summary of present lines
    # line_info = CSV.read(GRASS.datdir * "line_info.csv", DataFrame)

    # # get summary of present line atomic params etc.
    lp = LineProperties(exclude=[""])
    files = lp.file
    airwavs = GRASS.get_rest_wavelength(lp)
    names = GRASS.get_name(lp)

    # # number of breaks in plot needed
    # nregions = length(unique(files))
    # files = unique(files)
    # fdirs = split.(getindex.(splitdir.(files), 2), ".")

    # look for input data
    larsdir = "/storage/home/mlp95/ford_dir/mlp95/lars_spectra/"
    dirs = filter(isdir, readdir(larsdir, join=true))
    nregions = length(dirs)

    # get the wavelength regions
    wavmins = zeros(length(dirs))
    fluxmins = zeros(length(dirs))
    wavmaxs = zeros(length(dirs))
    fluxmaxs = zeros(length(dirs))
    for i in eachindex(dirs)
        # get subset'd dataframes
        df = GRASS.sort_spectrum_data(dir=dirs[i])
        df = subset(df, :axis => x -> x .== "c")
        file = joinpath(dirs[i], df.fname[1])

        # read in the file
        wavs, flux = GRASS.read_spectrum(file)
        flux = mean(flux,dims=2)./maximum(mean(flux,dims=2))

        wavmins[i] = minimum(wavs)
        fluxmins[i] = minimum(flux)
        wavmaxs[i] = maximum(wavs)
        fluxmaxs[i] = maximum(flux)
    end

    # sort on wavmin indices
    idx = sortperm(wavmins)
    wavmins .= wavmins[idx]
    fluxmins .= fluxmins[idx]
    wavmaxs .= wavmaxs[idx]
    fluxmaxs .= fluxmaxs[idx]
    dirs .= dirs[idx]

    # get width of wavelength regions
    wav_wids = wavmaxs .- wavmins
    wav_rats = wav_wids ./ wav_wids[1]

    # create plot objects
    fig = plt.figure(figsize=(21,8))
    gss = mpl.gridspec.GridSpec(1, nregions, width_ratios=wav_rats)
    axs = [plt.subplot(ax) for ax in gss]
    # fig, axs = plt.subplots(1, nregions, sharey=true, figsize=(15,3))
    fig.subplots_adjust(wspace=0.05)

    # loop over files
    for i in eachindex(dirs)
        # get subset'd dataframes
        df = GRASS.sort_spectrum_data(dir=dirs[i])
        df = subset(df, :axis => x -> x .== "c")
        file = joinpath(dirs[i], df.fname[1])

        # read in the file
        wavs, spec = GRASS.read_spectrum(file)

        # take a simple mean and roughly normalize
        wavs = dropdims(mean(wavs, dims=2), dims=2)
        spec = dropdims(mean(spec, dims=2), dims=2)
        spec ./= maximum(spec)

        # smooth the spectrum
        wavs = GRASS.moving_average(wavs, 5)
        spec = GRASS.moving_average(spec, 5)

        # interpolate up to high res for nice plotting of minima
        itp = GRASS.linear_interp(wavs, spec)
        wavs2 = range(minimum(wavs), maximum(wavs), step=minimum(diff(wavs)/5))
        spec2 = itp.(wavs2)

        # plot the data
        axs[i].plot(wavs2, spec2, c="k")

        # set the xlimits
        axs[i].set_xlim(minimum(wavs2) - 0.5, maximum(wavs2) + 0.5)
        axs[i].set_ylim(-0.1, 1.1)
        # axs[i].set_box_aspect(0.5)

        # axs[i].grid(false)

        # set the ticks
        wavmin = round(Int,minimum(wavs2))
        wavmid = round(Int, median(wavs2))
        wavmax = round(Int,maximum(wavs2))

        # find airwavs in wavelength region
        # idk = findall(x -> (x .<= maximum(wavs2) .& (x .>= minimum(wavs2))), airwavs)
        airwav_idx = findall(x -> (x .<= maximum(wavs2)) .& (x .>= minimum(wavs2)), airwavs)
        airwav_ann = airwavs[airwav_idx]
        names_ann = names[airwav_idx]

        texts = []
        for j in eachindex(airwav_idx)
            # find location on axis
            idx2 = findfirst(x-> x .>= airwav_ann[j], wavs2)
            min = argmin(spec2[idx2-50:idx2+50]) + idx2 - 50

            # set rotation
            if isapprox(airwav_ann[j], 5896, atol=1e0)
                rotation = 0.0
                d1 = 0.025
                d2 = 0.1
            else
                rotation = 270.0
                d1 = 0.025
                d2 = 0.22
            end

            # annotate with line name
            arrowprops = Dict("facecolor"=>"black", "lw"=>1.5, "arrowstyle"=>"-")
            txt = axs[i].annotate(("\${\\rm " * replace(names_ann[j], "_" => "\\ ") * "}\$"), rotation=rotation,
                                  (wavs2[min], spec2[min] - d1), (wavs2[min], spec2[min] - d2),
                                  arrowprops=arrowprops, horizontalalignment="center", fontsize=12)
            push!(texts, txt)
        end

        # d = 0.5
        # adjust_text.adjust_text(texts, arrowprops=Dict("arrowstyle"=>"-", "color"=>"k", "lw"=>0.5),
        #                         expand_text=(d, d), expand_points=(d, d), expand_objects=(d, d),
        #                         expand_align=(d,d))

        # set the xticks
        axs[i].set_xticks([wavmid-2, wavmid, wavmid+2])
        axs[i].xaxis.set_tick_params(rotation=45, fontsize=16)
        axs[i].yaxis.set_tick_params(fontsize=16)

        # deal with axis break decoration stuff
        d = 0.01
        if i > 1
            axs[i].yaxis.set_tick_params(left=false, labelleft=false)
            axs[i].spines["left"].set_visible(false)
            axs[i].plot([-d, +d], [-d, +d], transform=axs[i].transAxes, c="k", clip_on=false, lw=1)
            axs[i].plot([-d, +d], [1-d, 1+d], transform=axs[i].transAxes, c="k", clip_on=false, lw=1)
        end

        if i < length(dirs)
            axs[i].spines["right"].set_visible(false)
            axs[i].plot([1-d, 1+d], [-d, +d], transform=axs[i].transAxes, c="k", clip_on=false, lw=1)
            axs[i].plot([1-d, 1+d], [1-d, 1+d], transform=axs[i].transAxes, c="k", clip_on=false, lw=1)
        end
    end

    # make axis labels
    fig.supxlabel(L"{\rm Wavelength\ (\AA)}", y=-0.01, fontsize=21)
    fig.supylabel(L"{\rm Normalized\ Flux}", x=0.09, fontsize=21)

    # save the fig
    fig.savefig("spectra_collage.pdf")
    plt.clf(); plt.close()

    # plot all the lines with broken axes

    # # loop over lines in line list
    # for name in line_info.name
    #     # plot the spectrum
    #     plot_input_spectra(name, line_info; highlight=true)

    #     # plot the input data if it exists
    #     indir = GRASS.soldir * name * "/"
    #     if isdir(indir) && !isempty(indir)
    #         plot_input_data(name, line_info)
    #     end
    # end
# end

# if (run | plot)
#     main()
# end

