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

    # interpolate onto original lambda grid, extrapolate to continuum
    itp1 = linear_interp(allwavs, allints, bc=one(T))
    prof .*= itp1.(lambdas)
    return nothing
end

# fix for automated testing
line_profile! = line_profile_cpu!

# function to calc intensity at given x,y coord.
function line_profile_gpu!(mid::T, lambdas::AA{T,1}, prof::AA{T,1}, wsp::SynthWorkspace{T}) where T<:AF
    # move data to GPUU
    prof_gpu = CuArray(prof)
    lambdas_gpu = CuArray(prof)
    wavt_gpu = CuArray(wsp.wavt)
    dept_gpu = CuArray(wsp.dept)
    widt_gpu = CuArray(wsp.widt)
    lwavgrid_gpu = CuArray(lwavgrid)
    rwavgrid_gpu = CuArray(rwavgrid)
    allwavs_gpu = CuArray(allwavs)
    allints_gpu = CuArray(allints)

    # synthesize the line profile given bisector and width input data
    CUDA.@sync @cuda line_profile_gpu!(mid, lambdas_gpu, prof_gpu, wavt_gpu,
                                       dept_gpu, widt_gpu, lwavgrid_gpu,
                                       rwavgrid_gpu, allwavs_gpu, allints_gpu)

    prof .= Array(prof_gpu)
    return nothing
end

function line_profile_gpu!(mid, lambdas, prof, wavm, depm, widm, lwavgrid, rwavgrid, allwavs, allints)
    # set wavgrids to line center to start
    for i in 1:CUDA.length(lwavgrid)
        lwavgrid[i] = (mid - (0.5 * widm[i] - wavm[i]))
        rwavgrid[i] = (mid + (0.5 * widm[i] + wavm[i]))
    end
    rwavgrid[1] = lwavgrid[1] + 1e-3            # TODO: fix to deal with nodes

    # concatenate into one big array
    len = CUDA.length(rwavgrid)
    for i in 1:CUDA.length(rwavgrid)
        allwavs[i+len] = rwavgrid[i]
        allints[i+len] = depm[i]
        allwavs[i] = lwavgrid[CUDA.length(rwavgrid) - (i - 1)]
        allints[i] = depm[CUDA.length(rwavgrid) - (i - 1)]
    end

    # interpolate onto original lambda grid, extrapolate to continuum
    linear_interp_mult_gpu(prof, lambdas, allwavs, allints, 1.0)
    return nothing
end
