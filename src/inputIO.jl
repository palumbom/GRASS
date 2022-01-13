function sort_input_data(;dir::String=soldir, write::Bool=false)
    @assert isdir(dir)

    # get list of dirs containing data for each line
    linedirs = glob("*/", dir)

    # extract file parametersl put in data frame
    df = DataFrame(fpath=String[], fname=String[],
                   datetime=DateTime[], wave=String[],
                   mu=String[], axis=String[], ion=String[])

    # loop over directories
    for d in linedirs
        files = glob("*input.h5", d)
        if isempty(files)
            println(">>> No input data files found in " * d)
            continue
        end

        for f in files
            push!(df, extract_input_params(f))
        end
    end

    # sort dataframe in place
    sort!(df, [:axis, :mu, :datetime], rev=[false, true, false])

    # write to CSV
    if write
        CSV.write(df.wave[1] * "_bisectors.csv", df)
    end
    return df
end

function extract_input_params(s::String)
    # parse the fielname string
    paths = split(s, "/")
    fpath = join(paths[1:end-1], "/") * "/"
    fname = paths[end]
    fcomp = split(fname, "_")

    # parse out parameters
    species = fcomp[1]
    wavelength = fcomp[2]
    datetimes = DateTime(fcomp[3])
    muposition = fcomp[4]
    axis = fcomp[5]
    return [fpath fname datetimes wavelength muposition axis species]
end

function read_input_data(filename::String; masknans::Bool=false)
    wav = Array{Float64,2}[]
    bis = Array{Float64,2}[]
    wid = Array{Float64,2}[]
    dep = Array{Float64,2}[]
    h5open(filename, "r") do f
        g = f["input_data"]
        wav = read(g["wavelengths"])
        bis = read(g["bisectors"])
        wid = read(g["widths"])
        dep = read(g["depths"])
    end
    return wav, bis, wid, dep
end

function get_extension_dims(filename::String)
    out = 0
    h5open(filename, "r") do f
        g = f["input_data"]
        out = size(g["wavelengths"])
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
    widall = zeros(maximum(Nwav), sum(Ntim))
    depall = zeros(maximum(Nwav), sum(Ntim))

    # loop over the files, filling arrays
    for i in 1:Nfil
        wav, bis, wid, dep = read_input_data(df.fpath[i] * df.fname[i])
        wavall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= wav
        bisall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= bis
        widall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= wid
        depall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= dep
    end

    # adjust the mean if adjust_mean == true
    if adjust_mean
        wavall = adjust_data_mean(wavall, Ntim, Nfil)
        widall = adjust_data_mean(widall, Ntim, Nfil)
    end
    return wavall, bisall, widall, depall
end

function adjust_data_mean(array::AA{T,2}, Ntim::Vector{Int64}, Nfil::Int) where T<:Real
    # find the mean for the first chunk of bisectors
    group1 = array[:, 1:Ntim[1]]
    meangroup1 = mean(group1, dims=2)[:,1]
    for i in 2:Nfil
        # find the mean for the nth group of bisectors
        groupn = array[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])]
        meangroupn = mean(groupn, dims=2)[:,1]

        # find the distance between the means and correct by it
        meandist = meangroupn - meangroup1
        array[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .-= meandist
    end
    return array
end
