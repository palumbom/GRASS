function sort_spectrum_data(;dir::String="", write::Bool=false)
    # glob the files
    @assert isdir(dir);

    if isdir(dir * "spectra/")
        files = glob("*.fits", dir * "spectra/")
    else
        files = glob("*.fits", dir)
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

function write_line_params(line_df::DataFrame; clobber::Bool=false)
    # get the filename
    fname = GRASS.soldir * line_df.name[1] * ".h5"

    # create the file if it doesn't exist
    if clobber | !isfile(fname)
        h5open(fname, "w") do f; end
    end

    # write the line properties as attributes of the file
    h5open(fname, "r+") do f
        # get the attrributes
        attr = HDF5.attributes(f)

        # check if the metadata already exists
        if haskey(attr, "depth")
            println("\t >>> " * splitdir(fname)[end] * " metadata already exists...")
        else
            # read in the IAG data and isolate the line
            # TODO: line of interest may not be deepest line!!!
            println("\t >>> Writing line properties to " * splitdir(fname)[end])
            iag_wavs, iag_flux = read_iag_atlas(isolate=true, airwav=line_df.air_wavelength[1])

            idx1 = findfirst(x -> x .<= line_df.air_wavelength[1] - 0.25, iag_wavs)
            idx2 = findfirst(x -> x .>= line_df.air_wavelength[1] + 0.25, iag_wavs)
            iag_depth = 1.0 - minimum(view(iag_flux, idx1:idx2))

            # write the attributes to file metadata
            for n in names(line_df)
                if ismissing(line_df[!, n][1])
                    attr[n] = NaN
                else
                    attr[n] = line_df[!, n][1]
                end
            end
            attr["depth"] = iag_depth
        end

        # look for any pre-existing input data and delete it
        if !isempty(keys(f))
            println("\t >>> Purging old input data...")
            for k in keys(f)
                delete_object(f, k)
            end
        end
    end
    return nothing
end

function write_input_data(line_df::DataFrame, ax::String, mu::String,
                          datetime::Dates.DateTime, top_ints::AA{T,1},
                          bis::AA{T,2}, int1::AA{T,2}, wid::AA{T,2},
                          int2::AA{T,2}) where T<:AF
    # get the filename
    fname = GRASS.soldir * line_df.name[1] * ".h5"

    # create the file if it doesn't exist
    if !isfile(fname)
        h5open(fname, "w") do f; end
    end

    # write the data
    h5open(fname, "r+") do f
        # create the group for this disk position if it doesn't exists
        if !haskey(f, ax * "_" * mu)
            create_group(f, ax * "_" * mu)
            pos_group = f[ax * "_" * mu]

            # set attributes for this group
            attr = HDF5.attributes(pos_group)
            attr["mu"] = mu
            attr["axis"] = ax
        end

        # create the sub-group for this specific datetime
        pos_group = f[ax * "_" * mu]
        create_group(pos_group, string(datetime))
        g = pos_group[string(datetime)]

        # set attributes
        attr = HDF5.attributes(g)
        attr["datetime"] = string(datetime)
        attr["length"] = size(bis,2)

        # fill out the datasets
        g["top_ints"] = top_ints
        g["bisectors"] = bis
        g["bis_intensities"] = int1
        g["widths"] = wid
        g["wid_intensities"] = int2
    end
    return nothing
end

"""
BELOW IS LEGACY CODE
"""


function sort_bisector_data(;dir::String=soldir, write::Bool=false)
    # glob the files
    @assert isdir(dir)
    files = glob("*.bisect.fits", dir * "bisectors")
    @assert !isempty(files)

    # extract file parametersl put in data frame
    df = DataFrame(fpath=String[], fname=String[], datetime=DateTime[],
                   wave=String[], mu=String[], axis=String[], ion=String[])
    for i in 1:length(files)
        push!(df, extract_bisector_params(files[i]))
    end

    # sort dataframe in place
    sort!(df, [:axis, :mu, :datetime], rev=[false, true, false])

    # write to CSV
    if write
        CSV.write(df.wave[1] * "_bisectors.csv", df)
    end
    return df
end

# extract observation parameters from string name
function extract_bisector_params(s::String)
    # parse out filename
    paths = split(s, "/")
    fpath = join(paths[1:end-1], "/") * "/"
    fname = paths[end]
    fcomp = split(fname, "_")

    # parse out parameters
    datetimes = DateTime(fcomp[3], "YYYYmmdd-HHMMSS")
    wavelength = split(fcomp[4],"clv")[2]
    muposition = split(fcomp[5], ".")[1]
    if length(fcomp) == 7
        helio_axis = split(fcomp[6], ".")[1]
        ionspecies = split(fcomp[6], ".")[4]
    else
        helio_axis = "c"
        ionspecies = split(fcomp[5], ".")[4]
    end
    return [fpath fname datetimes wavelength muposition helio_axis ionspecies]
end

# simple function to read out wavelength/bisector vectors
function read_bisector(filename::String; masknans::Bool=false)
    wav = Array{Float64,1}[]
    bis = Array{Float64,1}[]
    FITS(filename) do f
        # λs in 1st extension, bis value in 2nd
        wav = FITSIO.read(f[1], :, :, 1)
        bis = FITSIO.read(f[1], :, :, 2)
    end

    # mask NaNs
    if masknans
        mask = .!isnan.(wav)
        wav = wav[mask]
        bis = bis[mask]
    end
    return wav, bis
end

function get_extension_dims(filename::String)
    out = 0
    FITS(filename) do f
        out = size(f[1])
    end
    return out
end

function get_number_waves(filename::String)
    dims = get_extension_dims(filename)
    return dims[1]
end

function get_number_times(filename::String)
    dims = get_extension_dims(filename)
    return dims[2]
end

# make one large time series for a given solar position
function stitch_time_series(df::DataFrame; adjust_mean::Bool=false, contiguous_only::Bool=false)
    # find out size of data
    if contiguous_only
        Nfil=1
    else
        Nfil = size(df, 1)
    end
    Ntim = map(get_number_times, df.fpath[i] * df.fname[i] for i in 1:Nfil)
    Nwav = map(get_number_waves, df.fpath[i] * df.fname[i] for i in 1:Nfil)

    # allocate array
    wavall = zeros(maximum(Nwav), sum(Ntim))
    bisall = zeros(maximum(Nwav), sum(Ntim))

    # loop over the files, filling arrays
    for i in 1:Nfil
        wav, bis = read_bisector(df.fpath[i] * df.fname[i])
        wavall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= wav
        bisall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= bis
    end

    # adjust the mean if adjust_mean == true
    if adjust_mean
        # find the mean for the first chunk of bisectors
        wavgroup1 = wavall[:, 1:Ntim[1]]
        meanwav1 = mean(wavgroup1, dims=2)[:,1]
        for i in 2:Nfil
            # find the mean for the nth group of bisectors
            wavgroupn = wavall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])]
            meanwavn = mean(wavgroupn, dims=2)[:,1]

            # find the distance between the means and correct by it
            meandist = meanwavn - meanwav1
            wavall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .-= meandist
        end
    end
    return wavall, bisall
end
