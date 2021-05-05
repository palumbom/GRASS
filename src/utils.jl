"""
Author: Michael Palumbo & Others
Contact: mlp95@psu.edu
Purpose: Various ancillary functions.
"""

function ismononotonic(A::AA{T,1}) where T<:AF
    return (all(diff(A) .>= zero(T)) | all(diff(A) .<= zero(T)))
end

function strip_columns(A::AA{T,2}, cols::AA{Bool,1}) where T<:AF
    @assert length(cols) == size(A,2)
    return A[:, .!cols]
end

"""
https://stackoverflow.com/questions/50899973/indices-of-unique-elements-of-vector-in-julia
"""
function uniqueidx(x::AA{T}) where T
    uniqueset = Set{T}()
    ex = eachindex(x)
    idxs = Vector{eltype(ex)}()
    for i in ex
        xi = x[i]
        if !(xi in uniqueset)
            push!(idxs, i)
            push!(uniqueset, xi)
        end
    end
    return idxs
end

function isunique(x::AA{T}) where T
    uniqueset = Set{T}()
    ex = eachindex(x)
    idxs = Vector{Bool}()
    for i in ex
        xi = x[i]
        if !(xi in uniqueset)
            push!(idxs, true)
            push!(uniqueset, xi)
        else
            push!(idxs, false)
        end
    end
    return idxs
end

"""
Credit: traktofon @ https://discourse.julialang.org/t/findnearest-function/4143/4
"""
function searchsortednearest(a::AbstractVector{T}, x::T) where T
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
        return idx
    else
        return idx-1
    end
end

function searchsortednearest(x::T, a::AbstractVector{T}) where T
    return searchsortednearest(a, x)
end

"""
TODO: this docstring
calc rms about the observed mean
"""
function calc_rms(A::AA{T,1}) where T<:AF
    return sqrt(sum((A.-mean(A)).^2)/length(A))
end

function strip_nans_by_column(A::AA{T,1}) where T<:AF
    nans = isnan.(A)
    return A[.!vec(nans)]
end

function strip_nans_by_column(A::AA{T,2}) where T<:AF
    nancols = all(isnan.(A), dims=1)
    return A[:,.!vec(nancols)]
end

function strip_nans_by_row(A::AA{T,2}) where T<:AF
    nanrows = all(isnan.(A), dims=2)
    return A[.!vec(nanrows),:]
end

function findnearest(A::AA{T,1}, t) where T<:Real
    return findmin(abs.(A.-t))[2]
end
