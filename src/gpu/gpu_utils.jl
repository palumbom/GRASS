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
