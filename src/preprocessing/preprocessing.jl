# function to write line parameters file
function write_line_params(line_df::DataFrame; clobber::Bool=false)
    # get the filename
    line_dir = GRASS.soldir * line_df.name[1] * "/"
    prop_file =  line_dir * line_df.name[1] * "_line_properties.h5"

    # don't bother if the file already exists
    if isfile(prop_file) & !clobber
        println("\t >>> " * splitdir(prop_file)[end] * " already exists...")
        return nothing
    end

    # read in the IAG data and isolate the line
    iag_wavs, iag_flux = read_iag(isolate=true, airwav=line_df.air_wavelength[1])
    iag_depth = 1.0 - minimum(iag_flux)

    # write the line properties file
    println("\t >>> Writing " * line_df.name[1] * "_line_properties.h5")
    h5open(prop_file, "w") do fid
        create_group(fid, "properties")
        g = fid["properties"]

        # fill out attributes
        attr = HDF5.attributes(g)
        for n in names(line_df)
            if ismissing(line_df[!, n][1])
                attr[n] = NaN
            else
                attr[n] = line_df[!, n][1]
            end
        end
        attr["depth"] = iag_depth
    end
    return nothing
end

function write_input_data(line_name, air_wavelength, fparams, wav, bis, dep, wid)
    # create output file name
    new_file = line_name * "_" * string(fparams[3]) * "_" * fparams[5] * "_" * fparams[6] * "_input.h5"

    # write the input data to the file
    h5open(GRASS.soldir * line_name * "/" * new_file, "w") do fid
        # create the group
        create_group(fid, "input_data")
        g = fid["input_data"]

        # make attributes
        attr = HDF5.attributes(g)
        attr["datetime"] = string(fparams[3])
        attr["air_wavelength"] = air_wavelength

        # convert mu to number
        mu_num = []
        for ch in fparams[5]
            push!(mu_num, tryparse(Int64, string(ch)))
        end

        new_string = prod(string.(mu_num[.!isnothing.(mu_num)]))
        if new_string[1] == '1'
            attr["mu"] = 1.0
            attr["axis"] = "c"
        elseif new_string[1] == '0'
            attr["mu"] = parse(Float64, "0." * new_string[2:end])
            attr["axis"] = fparams[6]
        end

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


function fit_line_wings(wavs_iso, flux_iso)
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

function replace_line_wings(fit, wavst, fluxt, min, val; debug=false)
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
