struct DiskParams{T1<:AF, T2<:Integer}
    N::T2
    Nt::T2
    pole::Tuple{T1,T1,T1}
    u1::T1
    u2::T1
end

"""
    DiskParams(; N=256, Nt=50, pole=(0.0, 1.0, 0.0))

Construct a `DiskParams` composite type instance.

# Arguments
- `N::Integer=256`: Length of N*N spatial grid
- `Nt::Integer=50`: Number of 15-second snapshots
- `pole::Tuple{Float64, Float64, Float64}=(0.0, 1.0, 0.0)`: Unit vector specificying rotation axis direction. Default is equator-on.
"""
function DiskParams(;N=256, Nt=50, pole=(0.0, 1.0, 0.0), u1=0.4, u2=0.26)
    # ensure pole is unit vector
    if sum(pole) != one(pole[1])
        pole = pole./sqrt(sum(pole.^2))
    end

    return DiskParams(N, Nt, pole, u1, u2)
end
