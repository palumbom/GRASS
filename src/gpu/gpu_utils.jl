function searchsortednearest_gpu(a,x)
    idx = CUDA.searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>CUDA.length(a)); return CUDA.length(a); end
    if (a[idx]==x); return idx; end
    if ((a[idx]-x)^2 < (a[idx-1]-x)^2)
        return idx
    else
        return idx-1
    end
end
