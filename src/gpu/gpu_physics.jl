function calc_mu_gpu(x, y, z, O⃗)
    dp = x * O⃗[1] + y * O⃗[2] + z * O⃗[3]
    n1 = CUDA.sqrt(O⃗[1]^2.0 + O⃗[2]^2.0 + O⃗[3]^2.0)
    n2 = CUDA.sqrt(x^2.0 + y^2.0 + z^2.0)
    return dp / (n1 * n2)
end

function calc_mu_gpu(x, y, z, Ox, Oy, Oz)
    dp = x * Ox + y * Oy + z * Oz
    n1 = CUDA.sqrt(Ox^2.0 + Oy^2.0 + Oz^2.0)
    n2 = CUDA.sqrt(x^2.0 + y^2.0 + z^2.0)
    return dp / (n1 * n2)
end

function sphere_to_cart_gpu(ρ, ϕ, θ)
    # compute trig quantities
    sinϕ = CUDA.sin(ϕ)
    sinθ = CUDA.sin(θ)
    cosϕ = CUDA.cos(ϕ)
    cosθ = CUDA.cos(θ)

    # now get cartesian coords
    x = ρ * cosϕ * sinθ
    y = ρ * sinϕ
    z = ρ * cosϕ * cosθ
    return x, y, z
end

function sphere_to_cart_gpu_eclipse(ρ, ϕ, θ)
    # compute trig quantities
    sinϕ = CUDA.sin(ϕ)
    sinθ = CUDA.sin(θ)
    cosϕ = CUDA.cos(ϕ)
    cosθ = CUDA.cos(θ)

    # now get cartesian coords
    x = ρ * cosϕ * cosθ
    y = ρ * cosϕ * sinθ
    z = ρ * sinϕ
    return x, y, z
end

function rotate_vector_gpu(x0, y0, z0, R_x)
    # do dot product
    x1 = x0 * R_x[1,1] + y0 * R_x[1,2] + z0 * R_x[1,3]
    y1 = x0 * R_x[2,1] + y0 * R_x[2,2] + z0 * R_x[2,3]
    z1 = x0 * R_x[3,1] + y0 * R_x[3,2] + z0 * R_x[3,3]
    return x1, y1, z1
end

function rotation_period_gpu(ϕ, A, B, C)
    sinϕ = sin(ϕ)
    return 360.0/(A + B * sinϕ^2.0 + C * sinϕ^4.0) 
end

function calc_dA_gpu(ρs, ϕc, dϕ, dθ)
    return ρs^2.0 * CUDA.sin(π/2.0 - ϕc) * dϕ * dθ
end

function quad_limb_darkening_gpu(μ, u1, u2)
    return 1.0 - u1 * (1.0 - μ) - u2 * (1.0 - μ)^2.0
end

function quad_limb_darkening_gpu(μ, u1, u2, u3, u4)
    return 1.0 - u1 * (1.0 - μ^0.5) - u2 * (1.0 - μ) - u3 * (1.0 - μ^1.5) - u4 * (1.0 - μ^2.0)
end

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

#---------------- here down needs clearing
function lorentzian_phase_curve(theta)
    # B_CB = 0.45
    # h_CB = 0.00055

    # x = tan(theta / 2) / h_CB
    # delta_CB = (B_CB/2) * ((1 + ((1- exp(-x)) / (x))) / ((1 + x)^2)) 
    # return (1 + delta_CB)

    B_SH = 0.55
    h_SH = 0.00019
    delta_SH = B_SH / (1 + (1/h_SH) * tan(theta / 2)) 
    return (1 + delta_SH)

    # return ((1 + delta_SH) + (1 + delta_CB)) / 2
end