# function to calc intensity at given x,y coord.
function line_profile!(i::T, j::T, mid::T, blueshift::T, lambdas::AA{T,1},
                       prof::AA{T,1}, wsp::SynthWorkspace{T};
                       pole::Tuple{T,T,T}=(zero(T),one(T),zero(T))) where T<:AF
    # calculate line center given rot. and conv. doppler shift -> λrest * (1 + z)
    λΔD = mid * (one(T) + patch_velocity_los(i, j, pole=pole) + blueshift/c_ms)

    # synthesize the line profile given bisector and width input data
    line_from_bis!(λΔD, lambdas, prof, wsp.wavt, wsp.dept, wsp.widt,
                   wsp.lwavgrid, wsp.rwavgrid, wsp.allwavs, wsp.allints)
    return nothing
end

function line_from_bis!(mid::T, lambdas::AA{T,1}, prof::AA{T,1},
                        wavm::AA{T,1}, depm::AA{T,1}, widm::AA{T,1},
                        lwavgrid::AA{T,1}, rwavgrid::AA{T,1},
                        allwavs::AA{T,1}, allints::AA{T,1}) where T<:AF
    # set wavgrids to line center to start
    lwavgrid .= (mid .- (0.5 .* widm .- wavm))
    rwavgrid .= (mid .+ (0.5 .* widm .+ wavm) .- (lwavgrid[2] - lwavgrid[1]))

    # concatenate into one big array
    len = length(rwavgrid)
    itr = len:-1:1
    allwavs[(len+1):end] .= rwavgrid
    allints[(len+1):end] .= depm
    allwavs[1:len] .= view(lwavgrid, itr)
    allints[1:len] .= view(depm, itr)

    # interpolate onto original lambda grid, extrapolate to continuum
    xs = (allwavs,)
    ys = allints
    it = Gridded(Linear())
    itp1 = extrapolate(interpolate!(T, xs, ys, it), 1.0)
    prof .*= itp1.(lambdas)
    return nothing
end



