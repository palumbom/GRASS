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
                          top_ints::AA{T,1}, bis::AA{T,2}, int::AA{T,2}, wid::AA{T,2}) where T<:AF
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
        g["intensities"] = int
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


function fit_line_wings(wavs_iso::AA{T,1}, flux_iso::AA{T,1}; debug::Bool=false) where T<:AF
    # get indices and values for minimum, depth, and bottom
    min = argmin(flux_iso)
    bot = flux_iso[min]
    depth = 1.0 - bot

    # get wing indices for various percentage depths into line
    lidx50, ridx50 = find_wing_index(0.5 * depth + bot, flux_iso, min=min)
    lidx90, ridx90 = find_wing_index(0.9 * depth + bot, flux_iso, min=min)

    # isolate the line wings and mask area around line core for fitting
    Δbot = 4
    core = min-Δbot:min+Δbot
    lwing = lidx90:lidx50
    rwing = ridx50:ridx90

    # create arrays to fit on
    wavs_lfit = vcat(wavs_iso[lwing], wavs_iso[core])
    flux_lfit = vcat(flux_iso[lwing], flux_iso[core])
    wavs_rfit = vcat(wavs_iso[core], wavs_iso[rwing])
    flux_rfit = vcat(flux_iso[core], flux_iso[rwing])

    # set boundary conditions and initial guess
    # GOOD FOR FeI 5434 + others
    if isapprox(wavs_iso[argmin(flux_iso)], 5896, atol=1e0)
        lb = [0.5, wavs_iso[min], 0.01, 0.05]
        ub = [2.5, wavs_iso[min], 0.75, 0.75]
        p0 = [.97, wavs_iso[min], 0.05, 0.16]
    else
        lb = [0.0, wavs_iso[min]-minimum(diff(wavs_iso))/2.0, 1e-4, 1e-4]
        ub = [2.5, wavs_iso[min]+minimum(diff(wavs_iso))/2.0, 0.15, 0.15]
        p0 = [0.3, wavs_iso[min], 0.02, 0.05]
    end

    # perform the fit
    lfit = curve_fit(GRASS.fit_voigt, wavs_lfit, flux_lfit, p0, lower=lb, upper=ub)
    rfit = curve_fit(GRASS.fit_voigt, wavs_rfit, flux_rfit, p0, lower=lb, upper=ub)

    if debug
        @show lfit.param
        @show rfit.param
        @show lfit.converged
        @show rfit.converged
    end
    return lfit, rfit
end

function replace_line_wings(lfit, rfit, wavst::AA{T,1}, fluxt::AA{T,1}, min::Int, val::T; debug::Bool=false) where T<:AF
    # get line model for all wavelengths in original spectrum
    lflux_new = GRASS.fit_voigt(wavst, lfit.param)
    rflux_new = GRASS.fit_voigt(wavst, rfit.param)

    # do a quick "normalization"
    lflux_new ./= maximum(lflux_new)
    rflux_new ./= maximum(rflux_new)

    # find indices
    idxl, idxr = find_wing_index(val, fluxt, min=min)

    # replace wings with model
    fluxt[1:idxl] .= lflux_new[1:idxl]
    fluxt[idxr:end] .= rflux_new[idxr:end]
    return nothing
end
