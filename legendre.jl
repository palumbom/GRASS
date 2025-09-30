using LinearAlgebra

"""
    legendreP(l, x)

Compute Legendre polynomial P_l(x) using recurrence.
"""
function legendreP(l::Int, x::Float64)
    if l == 0
        return 1.0
    elseif l == 1
        return x
    else
        Pnm2 = 1.0
        Pnm1 = x
        Pn = 0.0
        for n in 2:l
            Pn = ((2n - 1) * x * Pnm1 - (n - 1) * Pnm2) / n
            Pnm2, Pnm1 = Pnm1, Pn
        end
        return Pn
    end
end

"""
    legendre_for_l(l, θ_deg)

Compute normalized Legendre polynomial and its derivative
for a single l and colatitude θ (degrees).
"""
function legendre_for_l(l::Int, θ_deg::Float64)
    θ = deg2rad(θ_deg)
    cost = cos(θ)
    sint = sin(θ)

    # geodesy normalization
    norm_sht = sqrt(2 * l + 1)

    # l(l+1) normalization
    norm_l = sqrt(l * (l + 1))
    norm_l = norm_l == 0 ? 1.0 : norm_l  # avoid div by 0 for l=0

    # evaluate P_l(cost)
    P = legendreP(l, cost)

    # derivative dP_l/dz
    if l == 0
        dP = 0.0
    else
        Plm1 = legendreP(l - 1, cost)
        dP = l * (Plm1 - cost * P) / (1 - cost^2)
    end

    # apply normalizations
    leg     = P * norm_sht
    leg_dz  = dP * norm_sht

    leg_out    = leg / sqrt(2) / norm_l
    leg_d1_out = -sint * leg_dz / sqrt(2) / norm_l

    return leg_out, leg_d1_out
end

"""
    meridional_terms(θ_deg, ϕ_deg, B0_deg)

Compute the meridional circulation basis functions (s=2,4)
at a single latitude θ, longitude ϕ, and solar tilt B0.
All angles in degrees.
"""
function meridional_terms(θ_deg::Float64, ϕ_deg::Float64, B0_deg::Float64)
    θ  = deg2rad(θ_deg)
    ϕ  = deg2rad(ϕ_deg)
    B0 = deg2rad(B0_deg)

    cosB0, sinB0 = cos(B0), sin(B0)
    cosθ,  sinθ  = cos(θ),  sin(θ)
    cosϕ,  sinϕ  = cos(ϕ),  sin(ϕ)

    # geometric projection factor
    lt = sinB0 * sinθ - cosB0 * cosθ * cosϕ

    # compute Legendre derivative only for l=2 and l=4
    _, dt_l2 = legendre_for_l(2, θ_deg)
    _, dt_l4 = legendre_for_l(4, θ_deg)

    term_s2 = dt_l2 * lt
    term_s4 = dt_l4 * lt

    return term_s2, term_s4
end

# Example usage
t2, t4 = meridional_terms(30.0, 45.0, 5.94)
println("Meridional terms: s=2 → $t2, s=4 → $t4")
