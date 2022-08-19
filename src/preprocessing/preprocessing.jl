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
            println("\t >>> Writing line properties to " * splitdir(fname)[end])
            iag_wavs, iag_flux = read_iag(isolate=true, airwav=line_df.air_wavelength[1])
            iag_depth = 1.0 - minimum(iag_flux)

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

function write_input_data(line_df::DataFrame, ax::String, mu::String, datetime::Dates.DateTime,
                          wav::AA{T,2}, bis::AA{T,2}, dep::AA{T,2}, wid::AA{T,2}) where T<:AF#; clobber::Bool=False)
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

            # create the attributes for position group
            attr = HDF5.attributes(pos_group)
            attr["mu"] = mu
            attr["axis"] = ax
        end

        # create the sub-group for this specific datetime
        pos_group = f[ax * "_" * mu]
        create_group(pos_group, string(datetime))
        g = pos_group[string(datetime)]

        # fill out the datasets
        g["wavelengths"] = wav
        g["bisectors"] = bis
        g["depths"] = dep
        g["widths"] = wid
    end
    return nothing
end

function find_wing_index(val, arr; min=argmin(arr))
    lidx = findfirst(x -> x .>= val, reverse(arr[1:min]))
    ridx = findfirst(x -> x .>= val, arr[min:end])
    if isnothing(lidx)
        lidx = lastindex(reverse(arr[1:min])) - 1
    elseif isnothing(ridx)
        ridx = lastindex(arr[min:end]) - 1
    end
    return clamp(min - lidx, 1, min), clamp(ridx + min, min, length(arr))
end


function fit_line_wings(wavs_iso::AA{T,1}, flux_iso::AA{T,1}) where T<:AF
    # get indices and values for minimum, depth, and bottom
    min = argmin(flux_iso)
    bot = flux_iso[min]
    depth = 1.0 - bot

    # get wing indices for various percentage depths into line
    lidx50, ridx50 = find_wing_index(0.5 * depth + bot, flux_iso, min=min)
    lidx60, ridx60 = find_wing_index(0.6 * depth + bot, flux_iso, min=min)
    lidx70, ridx70 = find_wing_index(0.7 * depth + bot, flux_iso, min=min)
    lidx80, ridx80 = find_wing_index(0.8 * depth + bot, flux_iso, min=min)
    lidx90, ridx90 = find_wing_index(0.9 * depth + bot, flux_iso, min=min)

    # isolate the line wings and mask area around line core for fitting
    Δbot = 2
    core = min-Δbot:min+Δbot
    lwing = lidx90:lidx50
    rwing = ridx50:ridx90
    wavs_fit = vcat(wavs_iso[lwing], wavs_iso[core], wavs_iso[rwing])
    flux_fit = vcat(flux_iso[lwing], flux_iso[core], flux_iso[rwing])

    # set boundary conditions and initial guess
    # GOOD FOR FeI 5434 + others
    lb = [0.0, wavs_iso[min], 0.0, 0.0]
    ub = [1.0, wavs_iso[min], 0.5, 0.5]
    p0 = [1.0 - depth, wavs_iso[min], 0.02, 0.01]
    # GOOD FOR FeI 5434 + others

    # perform the fit
    fit = curve_fit(GRASS.fit_voigt, wavs_fit, flux_fit, p0, lower=lb, upper=ub)
    return fit
end

function replace_line_wings(fit, wavst::AA{T,1}, fluxt::AA{T,1}, min::Int, val::T; debug::Bool=false) where T<:AF
    # get line model for all wavelengths in original spectrum
    flux_new = GRASS.fit_voigt(wavst, fit.param)

    # do a quick "normalization"
    flux_new ./= maximum(flux_new)

    # find indices
    idxl, idxr = find_wing_index(val, fluxt, min=min)

    # adjust any "kinks" between the wing model and and raw spec
    flux_new_l = copy(flux_new)
    Δfluxl = fluxt[idxl] - flux_new_l[idxl]
    while Δfluxl > 0.0 && !isapprox(flux_new_l[idxl], fluxt[idxl], atol=1e-4)
        flux_new_l = circshift(flux_new_l, 1)
        Δfluxl = fluxt[idxl] - flux_new_l[idxl]
    end

    flux_new_r = copy(flux_new)
    Δfluxr = flux_new_r[idxr] - fluxt[idxr]
    while Δfluxr > 0.0 && !isapprox(flux_new_r[idxr], fluxt[idxr], atol=1e-4)
        flux_new_r = circshift(flux_new_r, -1)
        Δfluxr = fluxt[idxr] - flux_new_r[idxr]
    end

    # replace wings with model
    fluxt[1:idxl] .= flux_new_l[1:idxl]
    fluxt[idxr:end] .= flux_new_r[idxr:end]
    return nothing
end
