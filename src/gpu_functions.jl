function linear_interp_gpu(out, new_xs, xs, ys, bc)
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

function linear_interp_mult_gpu(out, new_xs, xs, ys, bc)
    # get GPU dims
    ix = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sx = blockDim().x * gridDim().x

    # perform the interpolation
    # n = CUDA.length(new_xs)
    # for i in 1:CUDA.length(new_xs)
    for i in ix:sx:length(new_xs)
        if (((new_xs[i] < CUDA.first(xs)) | (new_xs[i] > CUDA.last(xs))) & !CUDA.isnan(bc))
            out[i] = out[i] * bc
        elseif new_xs[i] <= CUDA.first(xs)
            out[i] = out[i] * CUDA.first(ys)
        elseif new_xs[i] >= CUDA.last(xs)
            out[i] = out[i] * CUDA.last(ys)
        else
            j = CUDA.searchsortedfirst(xs, new_xs[i]) - 1
            j0 = CUDA.clamp(j, CUDA.firstindex(ys), CUDA.lastindex(ys))
            j1 = CUDA.clamp(j+1, CUDA.firstindex(ys), CUDA.lastindex(ys))
            out[i] = out[i] * (ys[j0] + (ys[j1] - ys[j0]) * (new_xs[i] - xs[j0]) / (xs[j1] - xs[j0]))
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

function line_loop_gpu(prof, lines, depths, rot_shift, conv_blueshifts, lambdas, wavt, bist, dept, widt, lwavgrid, rwavgrid, allwavs, allints, top)
    # get GPU dims
    ix = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sx = blockDim().x * gridDim().x

    # synthesize the line
    for i in ix:sx:CUDA.length(lines)
        # slice out the data
        wavt1 = CUDA.view(wavt, :, i)
        bist1 = CUDA.view(bist, :, i)
        dept1 = CUDA.view(dept, :, i)
        widt1 = CUDA.view(widt, :, i)

        # slice out the correct spot in memory
        lwavgrid1 = CUDA.view(lwavgrid, :, i)
        rwavgrid1 = CUDA.view(rwavgrid, :, i)
        allwavs1 = CUDA.view(allwavs, :, i)
        allints1 = CUDA.view(allints, :, i)

        # first trim the bisectors to the correct depth
        trim_bisector_chop_gpu!(depths[i], wavt1, bist1, dept1, widt1, top)

        # calculate line center given rot. and conv. doppler shift -> λrest * (1 + z)
        λΔD = lines[i] * (1.0 + rot_shift) * (1.0 + conv_blueshifts[i])

        # update the line profile in place
        line_profile_gpu!(λΔD, lambdas, prof, wavt1, dept1, widt1, lwavgrid1, rwavgrid1, allwavs1, allints1)
    end
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

# function time_loop_gpu(t_loop::Int, prof::AA{T,1}, rot_shift::T,
#                        key::Tuple{Symbol, Symbol}, liter::UnitRange{Int},
#                        spec::SpecParams{T}, wsp::SynthWorkspace{T}; top::T=NaN) where T<:AF
function time_loop_gpu(t_loop, prof, rot_shift, key, liter,spec, wsp; top=NaN)
    # some assertions
    # @assert all(prof .== one(T))

    # get views needed for line synthesis (+ move data to GPU)
    wsp.wavt .= CUDA.view(spec.soldata.wav[key], :, t_loop)
    wsp.bist .= CUDA.view(spec.soldata.bis[key], :, t_loop)
    wsp.dept .= CUDA.view(spec.soldata.dep[key], :, t_loop)
    wsp.widt .= CUDA.view(spec.soldata.wid[key], :, t_loop)

    # # TODO figure this out
    # for i in eachindex(spec.variability)
    #     wsp.wavt[:,i] .*= spec.variability[i]
    # end

    # send the job to the gpu
    prof .= 1.0
    @cuda line_loop_gpu(prof, spec.lines, spec.depths, rot_shift,
                        spec.conv_blueshifts, spec.lambdas, wsp.wavt,
                        wsp.bist, wsp.dept, wsp.widt, wsp.lwavgrid,
                        wsp.rwavgrid, wsp.allwavs, wsp.allints, top)
    return nothing
end

"""
AF = AbstractFloat
AA = AbstractArray
using Pkg; Pkg.activate(".")
using CUDA
using GRASS
using JLD2
using FileIO
using DataFrames
using BenchmarkTools

# get data
data = GRASS.SolarData()
wavt_main = data.wav[(:c, :mu10)][:,1]
bist_main = data.bis[(:c, :mu10)][:,1]
dept_main = data.dep[(:c, :mu10)][:,1]
widt_main = data.wid[(:c, :mu10)][:,1]

# values
top = NaN
dep = 0.75
rot_shift = 0.0
conv_blueshift = 0.0
mid = 5434.5

# CPU stuff
spec_cpu = SpecParams(lines=[5434.5], depths=[0.75])
wsp_cpu = GRASS.SynthWorkspace(spec_cpu);
wavt_cpu = copy(wavt_main); wsp_cpu.wavt .= wavt_cpu
bist_cpu = copy(bist_main); wsp_cpu.bist .= bist_cpu
dept_cpu = copy(dept_main); wsp_cpu.dept .= dept_cpu
widt_cpu = copy(widt_main); wsp_cpu.widt .= widt_cpu

lambdas_cpu = range(mid-2.0, mid+2.0, step=mid/7e5);
prof_cpu = ones(length(lambdas_cpu));
allwavs_cpu = zeros(200);
allints_cpu = zeros(200);
lwavgrid_cpu = zeros(100);
rvwavgrid_cpu = zeros(100);

function bmark_cpu(lambdas, prof, wsp, spec)
    GRASS.time_loop_cpu(1, prof, 0.0, (:c, :mu10), 1:1, spec, wsp; top=NaN)
    return nothing
end

# GPU stuff
spec_gpu = SpecParams(lines=[5434.5], depths=[0.75])
wsp_gpu = GRASS.SynthWorkspace(spec_gpu);

wavt_gpu = CuArray(repeat(wavt_main, 1, length(spec_gpu.lines)))
bist_gpu = CuArray(repeat(bist_main, 1, length(spec_gpu.lines)))
dept_gpu = CuArray(repeat(dept_main, 1, length(spec_gpu.lines)))
widt_gpu = CuArray(repeat(widt_main, 1, length(spec_gpu.lines)))

lambdas_gpu = CuArray(range(mid-2.0, mid+2.0, step=mid/7e5));
prof_gpu = CUDA.ones(Float64, length(lambdas_gpu));
allwavs_gpu = CUDA.zeros(Float64, 200, CUDA.length(spec_gpu.lines));
allints_gpu = CUDA.zeros(Float64, 200, CUDA.length(spec_gpu.lines));
lwavgrid_gpu = CUDA.zeros(Float64, 100, CUDA.length(spec_gpu.lines));
rwavgrid_gpu = CUDA.zeros(Float64, 100, CUDA.length(spec_gpu.lines));

function bmark_gpu(lambdas, prof, wsp, spec)
    CUDA.@sync begin
        time_loop_gpu(1, prof, 0.0, (:c, :mu10), 1:1, spec, wsp, top=NaN)
    end
    return nothing
end

@btime bmark_cpu(lambdas_cpu, prof_cpu, wsp_cpu, spec_cpu)
@btime bmark_gpu(lambdas_gpu, prof_gpu, wsp_gpu, spec_gpu)

# do the CPU
# GRASS.trim_bisector_chop!(dep, wavt_cpu, bist_cpu, dept_cpu, widt_cpu, top=top)
# GRASS.line_loop(prof_cpu, mid, dep, rot_shift, conv_blueshift, lambdas_cpu, wsp)

# do the GPU
# @cuda trim_bisector_chop_gpu!(dep, wavt_gpu, bist_gpu, dept_gpu, widt_gpu, top)
# @cuda line_loop_gpu(prof_gpu, lines, depths, rot_shift, conv_blueshift, lambdas_gpu, wavt_gpu, bist_gpu, dept_gpu, widt_gpu, lwavgrid_gpu, rwavgrid_gpu, allwavs_gpu, allints_gpu, top)

# function bmark_cpu()

#     return
# end

# function bmark_gpu()


#     return
# end

# function compare(gpu, cpu)
#     println(maximum(abs.(Array(gpu) .- cpu)))
# end

# compare(wavt_gpu, wsp.wavt)
# compare(bist_gpu, wsp.bist)
# compare(dept_gpu, wsp.dept)
# compare(widt_gpu, wsp.widt)
# compare(prof_gpu, prof_cpu)
"""
