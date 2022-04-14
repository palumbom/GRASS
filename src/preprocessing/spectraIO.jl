function sort_spectrum_data(;dir::String=soldir, write::Bool=false)
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
