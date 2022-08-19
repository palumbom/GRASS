"""
function sort_input_data(;dir::String=soldir, write::Bool=false)
    @assert isdir(dir)

    # create data frame to store extracted paramaters
    df = DataFrame(fpath=String[], fname=String[],
                   datetime=DateTime[], wave=String[],
                   mu=String[], axis=String[], ion=String[])

    # get list of dirs containing data for each line
    linedirs = glob("*/", dir)

    # loop over directories
    if !isempty(linedirs)
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
    else
        files = glob("*input.h5", dir)
        if isempty(files)
            println(">>> No input data files found in " * dir)
            return nothing
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
"""

function sort_input_data(fname)
    # loop over directories
    if !isempty(linedirs)
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
    else
        files = glob("*input.h5", dir)
        if isempty(files)
            println(">>> No input data files found in " * dir)
            return nothing
        end
        for f in files
                push!(df, extract_input_params(f))
        end
    end

    # sort dataframe in place
    sort!(df, [:axis, :mu, :datetime], rev=[false, true, false])
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
    wav, bis, dep, wid = h5open(filename, "r") do f
        g = f["input_data"]
        wav = read(g["wavelengths"])
        bis = read(g["bisectors"])
        dep = read(g["depths"])
        wid = read(g["widths"])
        return wav, bis, dep, wid
    end
    return wav, bis, dep, wid
end

function get_extension_dims(filename::String)
    dims = h5open(filename, "r") do f
        g = f["input_data"]
        return size(g["wavelengths"])
    end
    return dims
end

function get_number_waves(filename::String)
    dims = get_extension_dims(filename)
    return dims[1]
end

function get_number_times(filename::String)
    dims = get_extension_dims(filename)
    return dims[2]
end

function parse_mu_string(s::String)
    s = s[3:end]
    return tryparse(Float64, s[1] * "." * s[2:end])
end

function parse_mu_string(s::Symbol)
    return parse_mu_string(string(s))
end

function parse_ax_string(s::String)
    if s == "c"; return 0; end;
    if s == "n"; return 1; end;
    if s == "s"; return 2; end;
    if s == "e"; return 3; end;
    if s == "w"; return 4; end;
end

function parse_ax_string(s::Symbol)
    return parse_ax_string(string(s))
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
    depall = zeros(maximum(Nwav), sum(Ntim))
    widall = zeros(maximum(Nwav), sum(Ntim))

    # loop over the files, filling arrays
    for i in 1:Nfil
        wav, bis, dep, wid = read_input_data(df.fpath[i] * df.fname[i])
        wavall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= wav
        bisall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= bis
        depall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= dep
        widall[:, sum(Ntim[1:i-1])+1:sum(Ntim[1:i])] .= wid
    end

    # adjust the mean if adjust_mean == true
    if adjust_mean
        wavall = adjust_data_mean(wavall, Ntim, Nfil)
        # widall = adjust_data_mean(widall, Ntim, Nfil)
    end
    return wavall, bisall, depall, widall
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
