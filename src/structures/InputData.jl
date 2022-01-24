struct InputData
    dirs::Array{String,1}
    lineprops::Array{LineProperties,1}
end

function InputData(dir::String, lineprops::LineProperties)
    return InputData([dir], [lineprops])
end

function InputData(;dir::String=soldir, kwargs...)
    # parse out filename info for all the input data
    @assert isdir(dir)
    df = sort_input_data(dir=dir)
    unique_dirs = unique(df.fpath)

    # fill in line properties
    lp = Array{LineProperties,1}(undef, length(unique_dirs))

    for i in eachindex(unique_dirs)
        prop_file = glob("*_line_properties.h5", unique_dirs[i])
        lp[i] = LineProperties(prop_file[1])
    end
    return InputData(unique_dirs, lp)
end

function InputData(config::String ;dir::String=soldir, kwargs...)



    return
end
