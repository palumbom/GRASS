# pull out DataFrames w/ file names and obs parameters
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


function sort_spectrum_data(;dir::String=soldir, write::Bool=false)
    # glob the files
    @assert isdir(dir);
    @assert isdir(dir * "spectra/")
    files = glob("*.chvtt.fits", dir * "spectra/")
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

function sort_width_data(;dir::String=soldir, write::Bool=false)
    # glob the files
    @assert isdir(dir)
    files = glob("*widths.fits", dir * "widths/")
    @assert !isempty(files)

    # extract file parametersl put in data frame
    df = DataFrame(fpath=String[], fname=String[], datetime=DateTime[],
                   wave=String[], mu=String[], axis=String[])
    for i in 1:length(files)
        push!(df, extract_width_params(files[i]))
    end

    # sort dataframe in place
    sort!(df, [:axis, :mu, :datetime], rev=[false, true, false])

    # write to CSV
    if write
        CSV.write(df.wave[1] * "_widths.csv", df)
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

function extract_width_params(s::String)
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
        # ionspecies = split(fcomp[6], ".")[4]
    else
        helio_axis = "c"
        # ionspecies = split(fcomp[5], ".")[4]
    end
    return [fpath fname datetimes wavelength muposition helio_axis]
end

function bin_spectrum(wavs::AbstractArray{T,2}, spec::AbstractArray{T,2};
                      binsize::Integer=10) where T<:Real
    @assert size(wavs) == size(spec)
    @assert binsize >= 1
    @assert binsize < size(wavs,2)

    # get sizes, allocate memory
    nwave = size(wavs,1)
    ntime = size(wavs,2)
    nbins = floor(Int, ntime/binsize)
    wavsb = zeros(nwave, nbins)
    specb = zeros(nwave, nbins)

    # do the binning
    for i in 0:(nbins-1)
        wavsb[:,i+1] = mean(wavs[:, i*binsize+1:(i+1)*binsize], dims=2)
        specb[:,i+1] = mean(spec[:, i*binsize+1:(i+1)*binsize], dims=2)
    end
    return wavsb, specb
end

function read_spectrum(filename::String)
    # Primary HDU has spectra + noise in frame 1/2
    spec = []
    nois = []
    head = []
    wavs = []
    FITS(filename) do f
        # read contents
        spec = FITSIO.read(f[1])[:,:,1]
        nois = FITSIO.read(f[1])[:,:,2]
        head = read_header(f[1])
        wavs = FITSIO.read(f[2]) #.* 1.0e10  # convert to angstroms

        # determine if conversion to angstroms is necessary
        if all(wavs .< 1e3)
            wavs .*= 1.0e10
        end
    end
    return wavs, convert.(Float64, spec)
end

function write_the_fits(fname::String, xdat::AbstractArray{T,2}, ydat::AbstractArray{T,2}) where T<:Real
    FITS(fname, "w") do io
        write(io, cat(xdat, ydat, dims=3))
    end
    return nothing
end
