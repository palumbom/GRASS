function measure_continuum(flux::AA{T,1}) where T<:AF
    # get initial guess
    continuum = maximum(flux)

    # iteratively remove flux measurements >std from continuum far from
    flux_temp = copy(flux)
    resid = continuum
    while abs(resid) > 1.0
        # set resi
        resid = continuum

        # get standard deviation and clip fluxes
        std_flux = std(flux_temp)
        idx_flux = abs.(flux_temp .- continuum) .<= 0.5 * std_flux
        flux_temp = flux_temp[idx_flux]

        # get peak of histogram
        h1 = fit(Histogram, flux_temp, nbins=100)
        idx1 = argmax(h1.weights)
        continuum = h1.edges[1][idx1]

        # compute condition
        resid -= continuum
    end

    # check if we overshot
    # if abs(maximum(flux) - continuum) < std(flux_temp)
        # return maximum(flux)
    # else
        return continuum
    # end
    # return nothing
end

function find_wing_index(val::T, arr::AA{T,1}; min::Int=argmin(arr)) where T<:AF
    @assert !isnothing(min)
    lidx = findfirst(x -> x .>= val, reverse(arr[1:min]))
    ridx = findfirst(x -> x .>= val, arr[min:end])
    if isnothing(lidx)
        lidx = lastindex(reverse(arr[1:min])) - 1
    elseif isnothing(ridx)
        ridx = lastindex(arr[min:end]) - 1
    end
    return clamp(min - lidx, 1, min), clamp(ridx + min, min, length(arr))
end

function bracket_line(λrest::T, wavs::AA{T,1}, flux::AA{T,1}; level::T=0.95) where T<:AF
    minbuff = 25
    lidx = findfirst(x -> x .>= λrest, wavs)
    lmin = argmin(flux[lidx-minbuff:lidx+minbuff]) + lidx - (minbuff+1)
    lbot = flux[lmin]
    depth = 1.0 - lbot

    # bracket the line
    idx1, idx2 = GRASS.find_wing_index(level * depth + lbot, flux, min=lmin)

    # take views of spectrum 
    wavs_iso = view(wavs, idx1:idx2)
    flux_iso = view(flux, idx1:idx2)
    return wavs_iso, flux_iso
end

function bracket_line(λrest::T, wavs::AA{T,1}, flux::AA{T,1}, nois::AA{T,1}; level::T=0.95) where T<:AF
    minbuff = 25
    lidx = findfirst(x -> x .>= λrest, wavs)
    lmin = argmin(flux[lidx-minbuff:lidx+minbuff]) + lidx - (minbuff+1)
    lbot = flux[lmin]
    depth = 1.0 - lbot

    # bracket the line
    idx1, idx2 = GRASS.find_wing_index(level * depth + lbot, flux, min=lmin)

    # take views of spectrum 
    wavs_iso = view(wavs, idx1:idx2)
    flux_iso = view(flux, idx1:idx2)
    nois_iso = view(nois, idx1:idx2)
    return wavs_iso, flux_iso, nois_iso
end

function fit_line_wings(wavs_iso::AA{T,1}, flux_iso::AA{T,1};
                        nois_iso::AA{T,1}=zeros(length(flux_iso)),
                        debug::Bool=false) where T<:AF
    # get indices and values for minimum, depth, and bottom
    min = argmin(flux_iso)
    bot = flux_iso[min]
    depth = 1.0 - bot

    if all(iszero.(nois_iso))
        println("derp")
        wts = ones(length(flux_iso))
    else
        wts = 1.0 ./ (nois_iso .^ 2.0)
    end

    # do first-pass fit to get ballpark lorentzian + gaussian comps
    p0 = [0.2, wavs_iso[min], 0.05, 0.05]
    lb = [0.0, wavs_iso[min]- 5.0 * minimum(diff(wavs_iso)), 1e-5, 1e-5]
    ub = [2.5, wavs_iso[min]+ 5.0 * minimum(diff(wavs_iso)), 0.15, 0.15]
    fit1 = curve_fit(GRASS.fit_voigt, wavs_iso, flux_iso, wts, p0, lower=lb, upper=ub)

    # get wing indices for various percentage depths into line
    lidx10, ridx10 = find_wing_index(0.10 * depth + bot, flux_iso, min=min)
    lidx20, ridx20 = find_wing_index(0.20 * depth + bot, flux_iso, min=min)
    lidx30, ridx30 = find_wing_index(0.30 * depth + bot, flux_iso, min=min)
    lidx40, ridx40 = find_wing_index(0.40 * depth + bot, flux_iso, min=min)
    lidx50, ridx50 = find_wing_index(0.50 * depth + bot, flux_iso, min=min)
    lidx60, ridx60 = find_wing_index(0.60 * depth + bot, flux_iso, min=min)
    lidx65, ridx65 = find_wing_index(0.65 * depth + bot, flux_iso, min=min)
    lidx70, ridx70 = find_wing_index(0.70 * depth + bot, flux_iso, min=min)
    lidx80, ridx80 = find_wing_index(0.80 * depth + bot, flux_iso, min=min)
    lidx85, ridx85 = find_wing_index(0.85 * depth + bot, flux_iso, min=min)
    lidx90, ridx90 = find_wing_index(0.90 * depth + bot, flux_iso, min=min)

    # isolate the line wings and mask area around line core for fitting
    Δbot = 1
    core = min-Δbot:min+Δbot
    atol = 0.5
    if isapprox(wavs_iso[argmin(flux_iso)], 5379.6, atol=0.5)
        lidx95, ridx95 = find_wing_index(0.95 * depth + bot, flux_iso, min=min)
        lwing = lidx90:lidx20
        rwing = ridx20:ridx95
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5432.9, atol=0.5)
        lidx95, ridx95 = find_wing_index(0.95 * depth + bot, flux_iso, min=min)
        lwing = lidx95:lidx20
        rwing = ridx20:ridx95
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5434.52, atol=0.3)
        lidx65, ridx65 = find_wing_index(0.65 * depth + bot, flux_iso, min=min)
        lwing = lidx65:lidx10
        rwing = ridx10:ridx65
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5435.8, atol=0.3)
        lwing = lidx90:lidx20
        rwing = ridx20:ridx90
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5436.3, atol=0.5)
        lidx90, ridx90 = find_wing_index(0.90 * depth + bot, flux_iso, min=min)
        lidx95, ridx95 = find_wing_index(0.95 * depth + bot, flux_iso, min=min)
        lwing = lidx90:lidx20
        rwing = ridx20:ridx95
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5380.3, atol=0.5)
        lwing = lidx70:lidx20
        rwing = ridx20:ridx70
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5383.3, atol=0.5)
        lidx90, ridx90 = find_wing_index(0.90 * depth + bot, flux_iso, min=min)
        lwing = lidx90:lidx40
        rwing = ridx40:ridx90
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5578.7, atol=0.5)
        lwing = lidx90:lidx30
        rwing = ridx30:ridx90
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5896.0, atol=0.5)
        lwing = lidx70:lidx40
        rwing = ridx40:ridx70
    elseif isapprox(wavs_iso[argmin(flux_iso)], 6149.25, atol=0.5)
        lwing = lidx90:lidx20
        rwing = ridx20:length(wavs_iso)
    elseif isapprox(wavs_iso[argmin(flux_iso)], 6151.62, atol=0.5)
        lidx98, ridx98 = find_wing_index(0.98 * depth + bot, flux_iso, min=min)
        lwing = lidx98:lidx30
        rwing = ridx30:ridx98
    elseif isapprox(wavs_iso[argmin(flux_iso)], 6170.5, atol=0.5)
        lidx95, ridx95 = find_wing_index(0.95 * depth + bot, flux_iso, min=min)
        lwing = lidx95:lidx10
        rwing = ridx20:ridx95
    elseif isapprox(wavs_iso[argmin(flux_iso)], 6173.3, atol=0.5)
        lidx95, ridx95 = find_wing_index(0.95 * depth + bot, flux_iso, min=min)
        lwing = lidx95:lidx10
        rwing = ridx10:ridx95
    elseif isapprox(wavs_iso[argmin(flux_iso)], 6301.5, atol=0.25)
        lwing = lidx90:lidx20
        rwing = ridx20:ridx90
    elseif isapprox(wavs_iso[argmin(flux_iso)], 6302.5, atol=0.25)
        lwing = lidx90:lidx20
        rwing = ridx20:ridx80
    else
        lwing = lidx80:lidx30
        rwing = ridx30:ridx80
    end

    # find ones
    idx_contl = findall(isone, flux_iso[1:min])
    idx_contr = findall(isone, flux_iso[min+1:end]) .+ min

    # create arrays to fit on
    wavs_lfit = vcat(wavs_iso[idx_contl], wavs_iso[lwing], wavs_iso[core])
    flux_lfit = vcat(flux_iso[idx_contl], flux_iso[lwing], flux_iso[core])

    wavs_rfit = vcat(wavs_iso[core], wavs_iso[rwing], wavs_iso[idx_contr])
    flux_rfit = vcat(flux_iso[core], flux_iso[rwing], flux_iso[idx_contr])

    # # create ararys to fit on
    # # don't demand we fit the core, in case of weird LTE departures
    # wavs_lfit = vcat(wavs_iso[idx_contl], wavs_iso[lwing])
    # flux_lfit = vcat(flux_iso[idx_contl], flux_iso[lwing])

    # wavs_rfit = vcat(wavs_iso[rwing], wavs_iso[idx_contr])
    # flux_rfit = vcat(flux_iso[rwing], flux_iso[idx_contr])

    if all(iszero.(nois_iso))
        wts_lfit = ones(length(flux_lfit))
        wts_rfit = ones(length(flux_rfit))
    else
        wts = 1.0 ./ (nois_iso .^ 2.0)
        # wts_lfit = vcat(nois_iso[idx_contl], nois_iso[lwing])
        # wts_rfit = vcat(nois_iso[rwing], nois_iso[idx_contr])

        wts_lfit = vcat(nois_iso[idx_contl], nois_iso[lwing], nois_iso[core])
        wts_rfit = vcat(nois_iso[core], nois_iso[rwing], nois_iso[idx_contr])
    end

    # set boundary conditions and initial guess

    if isapprox(wavs_iso[argmin(flux_iso)], 5896.0, atol=1e0)
        lb = [0.5, wavs_iso[min], 0.01, 0.05]
        ub = [2.5, wavs_iso[min], 0.75, 0.75]
        p0 = [.97, wavs_iso[min], 0.05, 0.16]
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5432.546, atol=0.25)
        lb = [0.0, wavs_iso[min]- 5.0 * minimum(diff(wavs_iso)), 1e-5, 1e-5]
        ub = [2.5, wavs_iso[min]+ 5.0 * minimum(diff(wavs_iso)), 0.15, 0.15]
        p0 = [0.2, wavs_iso[min], 0.05, 0.05]
    elseif isapprox(wavs_iso[argmin(flux_iso)], 5383.368, atol=0.25)
        lb = [0.0, wavs_iso[min]- 5.0 * minimum(diff(wavs_iso)), 1e-5, 1e-5]
        ub = [2.5, wavs_iso[min]+ 5.0 * minimum(diff(wavs_iso)), 0.15, 0.15]
        p0 = [0.2, wavs_iso[min], 0.05, 0.05]
    else
        lb = [0.0, wavs_iso[min]- 5.0 * minimum(diff(wavs_iso)), 1e-5, 1e-5]
        ub = [2.5, wavs_iso[min]+ 5.0 * minimum(diff(wavs_iso)), 0.15, 0.15]
        p0 = fit1.param
    end

    # perform the fit
    lfit = curve_fit(GRASS.fit_voigt, wavs_lfit, flux_lfit, wts_lfit, p0, lower=lb, upper=ub)
    rfit = curve_fit(GRASS.fit_voigt, wavs_rfit, flux_rfit, wts_rfit, p0, lower=lb, upper=ub)

    if debug
        @show fit1.param
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
