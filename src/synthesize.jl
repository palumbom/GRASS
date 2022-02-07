# function to calc intensity at given x,y coord.
function line_profile_cpu!(mid::T, lambdas::AA{T,1}, prof::AA{T,1}, wsp::SynthWorkspace{T}) where T<:AF
    # synthesize the line profile given bisector and width input data
    line_profile_cpu!(mid, lambdas, prof, wsp.wavt, wsp.dept, wsp.widt, wsp.lwavgrid, wsp.rwavgrid, wsp.allwavs, wsp.allints)
    return nothing
end

function line_profile_cpu!(mid::T, lambdas::AA{T,1}, prof::AA{T,1},
                           wavm::AA{T,1}, depm::AA{T,1}, widm::AA{T,1},
                           lwavgrid::AA{T,1}, rwavgrid::AA{T,1},
                           allwavs::AA{T,1}, allints::AA{T,1}) where T<:AF
    # set wavgrids to line center to start
    lwavgrid .= (mid .- (0.5 .* widm .- wavm))
    rwavgrid .= (mid .+ (0.5 .* widm .+ wavm))
    rwavgrid[1] = lwavgrid[1] + 1e-3            # TODO: fix to deal with nodes

    # concatenate into one big array
    len = length(rwavgrid)
    itr = len:-1:1
    allwavs[(len+1):end] .= rwavgrid
    allints[(len+1):end] .= depm
    allwavs[1:len] .= view(lwavgrid, itr)
    allints[1:len] .= view(depm, itr)


    # plt.plot(allwavs, allints);

    # interpolate onto original lambda grid, extrapolate to continuum
    itp1 = linear_interp(allwavs, allints, bc=one(T))
    prof .*= itp1.(lambdas)
    # plt.plot(lambdas, prof); plt.show()
    return nothing
end
