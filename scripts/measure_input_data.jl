using Pkg; Pkg.activate(".")
using CSV
using Glob
using DataFrames
using Statistics
using GRASS

import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using LaTeXStrings
mpl.style.use(GRASS.moddir * "fig.mplstyle")

# set LARS spectra absolute dir and read line info file
const data_dir = "/storage/group/ebf11/default/mlp95/lars_spectra/"
line_info = CSV.read(GRASS.datdir * "line_info.csv", DataFrame,
                           types=[String, String, Float64, Float64,
                                  Float64, Float64, Float64, Float64,
                                  String, Float64, Float64])

for col in eachcol(line_info)
    replace!(col,missing => NaN)
end

line_info[!,"height"] = convert.(Float64, line_info[!,"height"])
line_info[!,"avg_temp_80"] = convert.(Float64, line_info[!,"avg_temp_80"])
line_info[!,"avg_temp_50"] = convert.(Float64, line_info[!,"avg_temp_50"])

# output directories for misc stuff
const grassdir, plotdir, datadir = GRASS.check_plot_dirs()
if !isdir(plotdir * "spectra_fits")
    mkdir(plotdir * "spectra_fits")
end

function preprocess_line(line_name::String; clobber::Bool=true, verbose::Bool=true, debug::Bool=false)
    # delete old preprocessed data
    if !debug && clobber
        files = Glob.glob(line_name * ".h5", GRASS.soldir)
        rm.(files)
    end

    # find row with line info and write the line_params file
    line_df = subset(line_info, :name => x -> x .== line_name)
    if !debug
        GRASS.write_line_params(line_df)
    end

    # get the file name for output h5 file
    fname = GRASS.soldir * line_df.name[1] * ".h5"

    # parse out parameters from file names
    spec_df = GRASS.sort_spectrum_data(dir=data_dir * line_df.spectra_dir[1] * "/")
    sort!(spec_df, :mu, rev=true)

    # loop over the files
    fits_files = spec_df.fpath .* spec_df.fname
    for i in eachindex(fits_files)
        # debugging block + filename printing
        # if debug && i > 1
            # break
            # nothing
        if debug && splitdir(fits_files[i])[end] != "lars_l12_20161017-163508_clv6302_mu07_e.ns.chvtt.fits"
            continue
        # if debug && !contains(fits_files[i], "mu03_n")
            # continue
        elseif verbose
            println("\t >>> " * splitdir(fits_files[i])[end])
        end

        # get parameters for spectrum
        datetime = spec_df.datetime[i]
        mu_string = spec_df.mu[i]
        ax_string = spec_df.axis[i]

        # get nobs per bin
        cad = GRASS.get_observing_cadence(fits_files[i])
        if verbose
            println("\t \t >>> cadence = "  * string(cad) * " sec")
        end

        if !iszero(15 % cad)
            println("\t \t >>> Cadence cannot be binned to 15 seconds, moving on...")
            continue
        else
            binsize = convert(Int, 15/cad)
        end

        # read in the spectrum, and bin it to 15 second obs
        wavs, flux, nois = GRASS.bin_spectrum_weighted(GRASS.read_spectrum(fits_files[i])..., binsize=binsize)

        # calculation the continuum level
        # continua = GRASS.measure_continuum.(eachcol(flux))
        # for i in eachindex(continua)
        #     flux[:,i] ./= continua[i]
        #     nois[:,i] ./= continua[i]
        # end

        # normalize by continuum
        nois ./= maximum(flux, dims=1)
        flux ./= maximum(flux, dims=1)

        # allocate memory for measuring input data
        nflux = 100
        bis = zeros(nflux, size(wavs,2))
        int1 = zeros(nflux, size(wavs,2))
        int2 = zeros(nflux, size(wavs,2))
        wid = zeros(nflux, size(wavs,2))
        top = zeros(size(wavs,2))

        # loop over epochs in spectrum file
        for t in 1:size(wavs, 2)
            # debugging block
            # if debug && t > 1
                # break
            # if debug && t != 11
            #     continue
            # end

            # get view of this time slice
            wavst = view(wavs, :, t)
            fluxt = view(flux, :, t)
            noist = view(nois, :, t)

            if debug
                fig, (ax1, ax2) = plt.subplots(1,2, figsize=(9,6))
                ax1.plot(wavst, fluxt, c="k", label="raw spec")
            end

            # refine the location of the minimum
            minbuff = 25
            idx = findfirst(x -> x .>= line_df.air_wavelength[1], wavst)
            min = argmin(fluxt[idx-minbuff:idx+minbuff]) + idx - (minbuff+1)
            bot = fluxt[min]
            depth = 1.0 - bot
            if isnan(bot)
                println("\t\t >>> Problem with t = " * string(t) * ", moving on...")
                continue
            end

            # find indices to isolate the line
            idx1, idx2 = GRASS.find_wing_index(0.9 * depth + bot, fluxt, min=min)
            if isnothing(idx1) | isnothing(idx2)
                println("\t\t >>> Problem with t = " * string(t) * ", moving on...")
                continue
            end

            # check that the indices dont take us into another line
            if line_name == "FeI_6302"
                wavbuff = 0.2
            else
                wavbuff = 0.3
            end
            if line_name != "NaI_5896" && wavst[min] - wavst[idx1] > wavbuff
                idx1 = findfirst(x -> x .> wavst[min] - wavbuff, wavst)
                idx1 = argmax(fluxt[idx1:min]) + idx1
            end
            if line_name != "NaI_5896" && wavst[idx2] - wavst[min] > wavbuff
                idx2 = findfirst(x -> x .> wavst[min] + wavbuff, wavst)
                idx2 = argmax(fluxt[min:idx2]) + min
            end

            # pad flux_iso with ones
            cont_idxl = findall(x -> (1.0 .- x) .< 0.001, fluxt[1:idx1])
            cont_idxr = findall(x -> (1.0 .- x) .< 0.001, fluxt[idx2+2:end]) .+ idx2

            flux_padl = ones(length(cont_idxl))
            flux_padr = ones(length(cont_idxr))

            # view of isolated line
            wavs_iso = vcat(wavst[cont_idxl], copy(view(wavst, idx1:idx2)), wavst[cont_idxr])
            flux_iso = vcat(flux_padl, copy(view(fluxt, idx1:idx2)), flux_padr)
            nois_iso = vcat(noist[cont_idxl], copy(view(noist, idx1:idx2)), noist[cont_idxr])

            # if the spectrum is bad and line isolation didn't work, move on
            if line_name != "NaI_5896" && wavst[idx2] - first(wavst[idx1]) > 0.75
                println("\t\t >>> Problem with t = " * string(t) * ", moving on...")
                continue
            elseif std(flux_iso) < 0.01
                println("\t\t >>> Problem with t = " * string(t) * ", moving on...")
                continue
            end

            # fit the line wings
            lfit, rfit = GRASS.fit_line_wings(wavs_iso, flux_iso, nois_iso=nois_iso, debug=debug)
            if !(lfit.converged & rfit.converged)
                println("\t\t >>> Fit did not converge for t = " * string(t) * ", moving on...")
                if !debug
                    continue
                end
            elseif (abs(lfit.param[2] - line_df.air_wavelength[1]) > 0.1) | (abs(rfit.param[2] - line_df.air_wavelength[1]) > 0.1)
                println("\t\t >>> Wrong line was fit for t = " * string(t) * ", moving on...")
                if !debug
                    continue
                end
            end

            # pad the spectrum with ones if line is close to edge of spectral region
            specbuff = 1.5
            if wavst[min] - first(wavst) < specbuff
                wavs_pad = range(first(wavst) - specbuff, first(wavst), step=minimum(diff(wavst)))
                wavs_meas = vcat(wavs_pad[1:end-1], wavst)
                flux_meas = vcat(ones(length(wavs_pad)-1), fluxt)
                min += (length(wavs_pad) - 1)
            elseif last(wavst) - wavst[min] < specbuff
                wavs_pad = range(last(wavst), last(wavst) + specbuff, step=minimum(diff(wavst)))
                wavs_meas = vcat(wavst, wavs_pad[2:end])
                flux_meas = vcat(fluxt, ones(length(wavs_pad)-1))
            else
                wavs_meas = wavst
                flux_meas = fluxt
            end

            # replace the line wings above % continuum
            if line_name == "NaI_5896"
                top[t] = 0.65 * depth + bot
            elseif line_name == "FeI_5434"
                top[t] = 0.65 * depth + bot
            elseif line_name == "FeI_6302"
                top[t] = 0.7 * depth + bot
            elseif line_name == "CI_5380"
                top[t] = 0.7 * depth + bot
            else
                top[t] = 0.8 * depth + bot
            end

            # debugging code block
            if debug
                idxl, idxr = GRASS.find_wing_index(top[t], flux_meas, min=min)
                ax1.axhline(flux_meas[idxl], c="k", ls="--", alpha=0.5)
                ax1.axhline(flux_meas[idxr], c="k", ls="--", alpha=0.5)
                ax1.plot(wavs_meas[1:min], GRASS.fit_voigt(wavs_meas[1:min], lfit.param), c="tab:purple", label="left model")
                ax1.plot(wavs_meas[min:end], GRASS.fit_voigt(wavs_meas[min:end], rfit.param), c="tab:pink", label="right model")
            end

            # replace the line wings above top[t]% continuum
            GRASS.replace_line_wings(lfit, rfit, wavs_meas, flux_meas, min, top[t], debug=debug)

            # measure the bisector and width function
            bis[:,t], int1[:,t] = GRASS.calc_bisector(wavs_meas, flux_meas, nflux=nflux, top=0.999)
            int2[:,t], wid[:,t] = GRASS.calc_width_function(wavs_meas, flux_meas, nflux=nflux, top=0.999)

            # make sure the intensities are the same
            @assert all(int1[:,t] .== int2[:,t])

            # debug plottings
            if debug
                ax1.plot(wavs_meas, flux_meas, c="tab:orange", ls="--", label="cleaned")
                ax1.plot(bis[:,t], int1[:,t], c="tab:blue", label="bisector")
                ax1.set_xlabel("Wavelength")
                ax1.set_ylabel("Normalized Intensity")
                ax1.legend(fontsize=10)

                ax2.plot(int2[:,t], wid[:,t], c="tab:blue")
                ax2.axvline(flux_meas[idxl], c="k", ls="--", alpha=0.5)
                ax2.axvline(flux_meas[idxr], c="k", ls="--", alpha=0.5)
                ax2.set_xlabel("Normalized Intensity")
                ax2.set_ylabel("Width across line")
                fig.suptitle("\${\\rm " * replace(line_df.name[1], "_" => "\\ ") * "}\$")
                fig.savefig(plotdir * "spectra_fits/" * line_df.name[1] * ".pdf")
                plt.show()
                plt.clf(); plt.close()
            end
        end

        # remove columns that we skipped measuring
        bad_cols = zeros(Bool, size(int1,2))
        for i in 1:size(int1,2)
            if all(iszero.(int1))
                bad_cols[i] = true
            end
        end
        bis = GRASS.strip_columns(bis, bad_cols)
        int = GRASS.strip_columns(int1, bad_cols)
        wid = GRASS.strip_columns(wid, bad_cols)

        # write input data to disk
        if !debug
            GRASS.write_input_data(line_df, ax_string, mu_string, datetime, top, bis, int, wid)
        end
    end

    # now go back through and compute convective blueshift for each disk pos
    println("\t >>> Measuring convective blueshifts...")
    if !debug
        GRASS.measure_convective_blueshifts(fname)
    end
    return nothing
end

function main()
    for name in line_info.name
        # skip the "hard" lines for now
        # (name in ["CI_5380", "FeI_5382"]) && continue
        name != "FeI_6302" && continue

        # print the line name and preprocess it
        println(">>> Processing " * name * "...")
        preprocess_line(name, debug=false)
    end
    return nothing
end

main()
