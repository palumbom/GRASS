function linear_interp_gpu(out, new_xs, xs, ys, bc)
    # get GPU dims
    i = threadIdx().x + blockDim().x * (blockIdx().x-1)

    # perform the interpolation
    n = CUDA.length(new_xs)
    for i in 1:CUDA.length(new_xs)
        if (((new_xs[i] < CUDA.first(xs)) | (new_xs[i] > CUDA.last(xs))) & !CUDA.isnan(bc))
            out[i] = bc
        elseif new_xs[i] <= CUDA.first(xs)
            out[i] = CUDA.first(ys)
        elseif new_xs[i] >= CUDA.last(xs)
            out[i] = CUDA.last(ys)
        else
            j = CUDA.searchsortedfirst(xs, new_xs[i]) - 1
            j0 = CUDA.clamp(j, CUDA.firstindex(ys), CUDA.lastindex(ys))
            j1 = CUDA.clamp(j+1, CUDA.firstindex(ys), CUDA.lastindex(ys))
            out[i] = ys[j0] + (ys[j1] - ys[j0]) * (new_xs[i] - xs[j0]) / (xs[j1] - xs[j0])
        end
    end
    return nothing
end

function trim_bisector_chop_gpu!(depth, wavt, bist, dept, widt, top)
    # replace spurious measurements at top of bisector
    ind1 = CUDA.searchsortedfirst(bist, 1.0 - depth)

    # TODO this will kill the code
    if !CUDA.isnan(top)
        ind2 = CUDA.searchsortedfirst(bist, top)
        wavt[ind2:end] .= wavt[ind2]
    end

    # get knots
    xs = CUDA.view(bist, ind1:CUDA.length(bist))
    ys = CUDA.view(wavt, ind1:CUDA.length(wavt))

    # get new grid of depths
    step = depth/(CUDA.length(dept) - 1)
    for i in 1:CUDA.length(dept)
        dept[i] = (1.0 - depth) + (i-1) * step
    end

    # do the interpolation, assign results to memory and return
    linear_interp_gpu(wavt, dept, xs, ys, NaN)

    # now assign bisector fluxes from dept
    for i in 1:CUDA.length(bist)
        bist[i] = dept[i]
    end
    return nothing
end

function line_loop_gpu(prof, mid, depth, rot_shift, conv_blueshift, lambdas, wavt, bist, dept, widt, top)
    # first trim the bisectors to the correct depth
    trim_bisector_chop_gpu!(depth, wavt, bist, dept, widt, top)

    # calculate line center given rot. and conv. doppler shift -> λrest * (1 + z)
    λΔD = mid * (1.0 + rot_shift) * (1.0 + conv_blueshift)

    # find window around shifted line
    lind = CUDA.findfirst(x -> x > λΔD - 0.5, lambdas)
    rind = CUDA.findfirst(x -> x > λΔD + 0.5, lambdas)



    return nothing
end

function line_profile_gpu!(mid, lambdas, prof, wavt, dept, widt, lwavgrid, rwavgrid, allwavs, allints)


    return
end

using Pkg; Pkg.activate(".")
using CUDA
using GRASS
using JLD2
using FileIO
using DataFrames

data = GRASS.SolarData()
wavt_main = data.wav[(:c, :mu10)][:,1]
bist_main = data.bis[(:c, :mu10)][:,1]
dept_main = data.dep[(:c, :mu10)][:,1]
widt_main = data.wid[(:c, :mu10)][:,1]

wavt1 = CuArray(wavt_main)
bist1 = CuArray(bist_main)
dept1 = CuArray(dept_main)
widt1 = CuArray(widt_main)

const top = NaN
const dep = 0.75
const rot_shift = 0.0
const conv_blueshift = 0.0
const mid = 5434.5
lambdas = CuArray(range(mid-0.75, mid+0.75, step=mid/7e5))
prof = CuArray(ones(length(lambdas)))

# @cuda trim_bisector_chop_gpu!(dep, wavt1, bist1, dept1, widt1, NaN)
@cuda line_loop_gpu(prof, mid, dep, rot_shift, conv_blueshift, lambdas, wavt1, bist1, dept1, widt1, top)
