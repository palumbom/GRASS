# function to calc intensity at given x,y coord.
function line_profile_cpu!(mid::T, dA::T, ld::T, ext::T, lambdas::AA{T,1}, prof::AA{T,1}, wsp::SynthWorkspaceEclipse{T}, ext_toggle::Bool) where T<:AF
    # synthesize the line profile given bisector and width input data
    line_profile_cpu!(mid, dA, ld, ext, lambdas, prof, wsp.bist, wsp.intt, wsp.widt,
                      wsp.lwavgrid, wsp.rwavgrid, wsp.allwavs, wsp.allints, ext_toggle)
    return nothing
end