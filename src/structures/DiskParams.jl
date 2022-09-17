struct DiskParams{T<:AF}
    N::Int
    Nt::Int
    pole::Tuple{T,T,T}
    u1::T
    u2::T
end

"""
    DiskParams(; N=132, Nt=50, pole=(0.0, 1.0, 0.0))

Construct a `DiskParams` composite type instance.

# Arguments
- `N::Integer=132`: Length of N*N spatial grid
- `Nt::Integer=50`: Number of 15-second snapshots
- `pole::Tuple{Float64, Float64, Float64}=(0.0, 1.0, 0.0)`: Unit vector specificying rotation axis direction. Default is equator-on.
"""
function DiskParams(;N=132, Nt=NaN, pole=(0.0, 1.0, 0.0), u1=0.4, u2=0.26)
    # assertions and warnings
    @assert !isnan(Nt)
    if N != 132
        @warn "N should be set to 132 for physical validity!"
    end

    # ensure pole is unit vector
    if sum(pole) != one(pole[1])
        @warn "Given pole vector is not a unit vector. Normalizing..."
        pole = pole./sqrt(sum(pole.^2))
    end
    return DiskParams(N, Nt, convert.(typeof(u1), pole), u1, u2)
end
