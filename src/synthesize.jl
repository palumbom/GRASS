# function to calc intensity at given x,y coord.
function line_profile_cpu!(mid::T, dA::T, ld::T, ext::T, contrast::T, lambdas::AA{T,1}, prof::AA{T,1}, wsp::SynthWorkspaceEclipse{T}) where T<:AF
    # synthesize the line profile given bisector and width input data
    line_profile_cpu!(mid, dA, ld, ext, contrast, lambdas, prof, wsp.bist, wsp.intt, wsp.widt,
                      wsp.lwavgrid, wsp.rwavgrid, wsp.allwavs, wsp.allints)
    return nothing
end

function line_profile_cpu!(mid::T, dA::T, ld::T, ext::T, contrast::T, lambdas::AA{T,1}, prof::AA{T,1},
                           bism::AA{T,1}, intm::AA{T,1}, widm::AA{T,1},
                           lwavgrid::AA{T,1}, rwavgrid::AA{T,1},
                           allwavs::AA{T,1}, allints::AA{T,1}) where T<:AF
    # set wavgrids to line center to start
    lwavgrid .= (mid .- (0.5 .* widm .- bism))
    rwavgrid .= (mid .+ (0.5 .* widm .+ bism))

    # concatenate into one big array
    len = length(rwavgrid)
    itr = len:-1:1
    allwavs[(len+1):end] .= rwavgrid
    allints[(len+1):end] .= intm
    allwavs[1:len] .= view(lwavgrid, itr)
    allints[1:len] .= view(intm, itr)

    # get indices for interpolation view
    lind = findlast(x -> x <= allwavs[1], lambdas)
    if isnothing(lind)
        lind = firstindex(lambdas)
    end

    rind = findfirst(x -> x >= allwavs[end], lambdas)
    if isnothing(rind)
        rind = lastindex(lambdas)
    end

    # get views
    lambda_window = view(lambdas, lind:rind)
    prof_window = view(prof, lind:rind)

    # interpolate onto original lambda grid, extrapolate to continuum
    itp1 = linear_interp(allwavs, allints, bc=one(T))
    prof_window .+= itp1.(lambda_window) .* dA .* ld .* contrast #.* ext

    # make sure other lambda values get weight
    if lind != firstindex(lambdas)
        prof[1:lind-1] .+= (dA .* ld .* contrast) #.* ext
    end

    if rind != lastindex(lambdas)
        prof[rind+1:end] .+= (dA .* ld .* contrast) #.* ext
    end

    return nothing
end

# function to calc intensity at given x,y coord.
function line_profile_cpu!(mid::T, weight::T, lambdas::AA{T,1}, prof::AA{T,1}, wsp::SynthWorkspace{T}) where T<:AF
    # synthesize the line profile given bisector and width input data
    line_profile_cpu!(mid, weight, lambdas, prof, wsp.bist, wsp.intt, wsp.widt,
                      wsp.lwavgrid, wsp.rwavgrid, wsp.allwavs, wsp.allints)
    return nothing
end

function line_profile_cpu!(mid::T, weight::T, lambdas::AA{T,1}, prof::AA{T,1},
                           bism::AA{T,1}, intm::AA{T,1}, widm::AA{T,1},
                           lwavgrid::AA{T,1}, rwavgrid::AA{T,1},
                           allwavs::AA{T,1}, allints::AA{T,1}) where T<:AF
    # set wavgrids to line center to start
    lwavgrid .= (mid .- (0.5 .* widm .- bism))
    rwavgrid .= (mid .+ (0.5 .* widm .+ bism))

    # concatenate into one big array
    len = length(rwavgrid)
    itr = len:-1:1
    allwavs[(len+1):end] .= rwavgrid
    allints[(len+1):end] .= intm
    allwavs[1:len] .= view(lwavgrid, itr)
    allints[1:len] .= view(intm, itr)

    # get indices for interpolation view
    lind = findlast(x -> x <= allwavs[1], lambdas)
    if isnothing(lind)
        lind = firstindex(lambdas)
    end

    rind = findfirst(x -> x >= allwavs[end], lambdas)
    if isnothing(rind)
        rind = lastindex(lambdas)
    end

    # get views
    lambda_window = view(lambdas, lind:rind)
    prof_window = view(prof, lind:rind)

    # interpolate onto original lambda grid, extrapolate to continuum
    itp1 = linear_interp(allwavs, allints, bc=one(T))
    prof_window .+= itp1.(lambda_window) .* weight

    # make sure other lambda values get weight
    if lind != firstindex(lambdas)
        prof[1:lind-1] .+= weight
    end

    if rind != lastindex(lambdas)
        prof[rind+1:end] .+= weight
    end

    return nothing
end