struct InputData
    soldata::Array{SolarData,1}
    lp::Array{LineProperties,1}
end

function InputData(;dir::String=soldir, kwargs...)
    # parse out filename info for all the input data
    @assert isdir(dir)
    idf = sort_input_data(dir=dir)
    unique_dirs = unique(idf.fpath)

    # fill in line properties
    lp = Array{LineProperties,1}(undef, length(unique_dirs))

    for i in eachindex(unique_dirs)
        prop_file = glob("*_line_properties.h5", unique_dirs[i])
        lp[i] = LineProperties(prop_file[1])
    end

    # get the rest wavelengths
    λrest = get_rest_wavelength.(lp)

    # allocate memory for data and loop through directory structure
    soldata = Array{SolarData,1}(undef, length(lp))
    for i in eachindex(lp)
        soldata[i] = SolarData(dir=unique_dirs[i], λrest=λrest[i])
    end
    return InputData(soldata, lp)
end
