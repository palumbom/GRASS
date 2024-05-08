# alias the CUDA.@sync macro (stolen from CUDA.jl source)
macro cusync(ex...)
    # destructure the `@sync` expression
    code = ex[end]
    kwargs = ex[1:end-1]

    # decode keyword arguments
    for kwarg in kwargs
        Meta.isexpr(kwarg, :(=)) || error("Invalid keyword argument $kwarg")
        key, val = kwarg.args
        if key == :blocking
            Base.depwarn("the blocking keyword to @sync has been deprecated", :sync)
        else
            error("Unknown keyword argument $kwarg")
        end
    end

    quote
        local ret = $(esc(code))
        CUDA.synchronize()
        ret
    end
end

function searchsortednearest_gpu(a,x)
    idx = CUDA.searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>CUDA.length(a)); return CUDA.length(a); end
    if (a[idx]==x); return idx; end
    if (CUDA.abs(a[idx]-x) < CUDA.abs(a[idx-1]-x))
        return idx
    else
        return idx-1
    end
end

function filter_array_gpu!(output::AA{T,1}, input::AA{T,1}, pred, n) where T
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    for i in idx:sdx:CUDA.size(input,1)
        if CUDA.isone(pred[i])
            n += 1
            output[n] = input[i]
        end
    end
    return nothing
end

function filter_array_gpu!(output::AA{T,2}, input::AA{T,2}, pred, n) where T
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = gridDim().x * blockDim().x

    for i in idx:sdx:CUDA.size(input,1)
        if CUDA.isone(pred[i])
            n += 1
            for j in 1:CUDA.size(input, 2)
                output[n, j] = input[i, j]
                continue
            end
        end
    end
    return nothing
end

