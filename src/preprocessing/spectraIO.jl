function sort_spectrum_data(;dir::String="", write::Bool=false)
    # glob the files
    @assert isdir(dir);

    if isdir(dir * "spectra/")
        files = glob("*.chvtt.fits", dir * "spectra/")
    else
        files = glob("*.chvtt.fits", dir)
    end
    @assert !isempty(files)

    # extract file parametersl put in data frame
    df = DataFrame(fpath=String[], fname=String[], datetime=DateTime[],
                   wave=String[], mu=String[], axis=String[])
    for i in 1:length(files)
        push!(df, extract_line_params(files[i]))
    end

    # sort dataframe in place
    sort!(df, [:axis, :mu, :datetime], rev=[false, true, false])

    # write to CSV
    if write
        CSV.write(df.wave[1] * "_lines.csv", df)
    end
    return df
end

function extract_line_params(s::String)
    # parse out filename
    paths = split(s, "/")
    fpath = join(paths[1:end-1], "/") * "/"
    fname = paths[end]
    fcomp = split(fname, "_")

    # parse out parameters
    datetimes = DateTime(fcomp[3], "YYYYmmdd-HHMMSS")
    wavelength = split(fcomp[4],"clv")[2]
    muposition = split(fcomp[5], ".")[1]
    if length(fcomp) == 6
        helio_axis = split(fcomp[6], ".")[1]
    else
        helio_axis = "c"
    end
    return [fpath fname datetimes wavelength muposition helio_axis]
end

function get_observing_cadence(filename::String)
    # Primary HDU has spectra + noise in frame 1/2
    cad = FITS(filename) do f
        # read contents
        head = read_header(f[1])
        return head["CADENCE"]
    end
    return cad
end


function read_spectrum(filename::String; omit_bad_rows::Bool=true)
    # Primary HDU has spectra + noise in frame 1/2
    head = []
    flux = []
    nois = []
    wavs = []
    FITS(filename) do f
        # read contents
        head = read_header(f[1])
        flux = convert.(Float64, FITSIO.read(f[1])[:,:,1])
        nois = convert.(Float64, FITSIO.read(f[1])[:,:,2])
        wavs = FITSIO.read(f[2])

        # get table data
        if omit_bad_rows
            # find rows with bad status
            df = DataFrame(f[3])
            cols = map(iszero, hcat(df.COMB_F0, df.COMB_FR, df.COMB_LOCKED, df.COMB_PLO_OK))
            status = .!map(all, eachrow(cols))

            # remove the bad columns
            wavs = strip_columns(wavs, status)
            flux = strip_columns(flux, status)
            nois = strip_columns(nois, status)
        end

        # determine if conversion to angstroms is necessary
        if all(wavs .< 1e3)
            wavs .*= 1.0e10
        end
    end
    return wavs, flux, nois
end

function bin_spectrum(wavs::AA{T,2}, flux::AA{T,2}; binsize::Integer=10) where T<:Real
    @assert size(wavs) == size(flux)
    @assert binsize >= 1
    @assert binsize < size(wavs,2)

    # get sizes
    nwave = size(wavs,1)
    ntime = size(wavs,2)
    nbins = floor(Int, ntime/binsize)

    # do the binning
    for i in 0:(nbins-1)
        # get new wavelength grid
        new_wavs = dropdims(mean(wavs[:, i*binsize+1:(i+1)*binsize], dims=2), dims=2)

        # do-flux conserving re-bin onto new wavelength grid
        for j in i*binsize+1:(i+1)*binsize
            flux[:, j] = rebin_spectrum(view(wavs, :, j), view(flux, :, j), new_wavs)
        end

        # assign new wavelength to memory and take mean flux
        wavs[:,i+1] = new_wavs
        flux[:,i+1] = mean(flux[:, i*binsize+1:(i+1)*binsize], dims=2)
    end
    return view(wavs, :, 1:nbins), view(flux, :, 1:nbins)
end

function bin_spectrum_weighted(wavs::AA{T,2}, flux::AA{T,2},
                               nois::AA{T,2}; binsize::Integer=10) where T<:Real
    @assert size(wavs) == size(flux)
    @assert binsize >= 1
    @assert binsize < size(wavs,2)

    # get sizes
    nwave = size(wavs,1)
    ntime = size(wavs,2)
    nbins = floor(Int, ntime/binsize)

    # do the binning
    for i in 0:(nbins-1)
        # get new wavelength grid
        new_wavs = dropdims(mean(wavs[:, i*binsize+1:(i+1)*binsize], dims=2), dims=2)

        # do-flux conserving re-bin onto new wavelength grid
        for j in i*binsize+1:(i+1)*binsize
            flux_new, nois_new = rebin_spectrum(view(wavs, :, j), view(flux, :, j),
                                                view(nois, :, j), new_wavs)

            flux[:, j] = flux_new
            nois[:, j] = nois_new
        end

        # assign new wavelength to memory
        wavs[:,i+1] = new_wavs

        # get weighted flux
        wghts = 1.0 ./ view(nois, :,i*binsize+1:(i+1)*binsize).^2.0
        flux[:,i+1] = sum(flux[:, i*binsize+1:(i+1)*binsize] .* wghts, dims=2)
        flux[:,i+1] ./= sum(wghts, dims=2)

        # get the new uncertainty on weighted average
        nois[:,i+1] = 1.0 ./ sqrt.(sum(wghts, dims=2))
    end
    return view(wavs, :, 1:nbins), view(flux, :, 1:nbins), view(nois, :, 1:nbins)
end
