struct DiskParams{T1<:AF, T2<:Integer}
    N::T2
    Nt::T2
    pole::Tuple{T1,T1,T1}
    u1::T1
    u2::T1
end

function DiskParams(;N=128, Nt=25, pole=(0.0, 1.0, 0.0), u1=0.4, u2=0.26)
    # ensure pole is unit vector
    if sum(pole) != one(pole[1])
        pole = pole./sqrt(sum(pole.^2))
    end

    return DiskParams(N, Nt, pole, u1, u2)
end
