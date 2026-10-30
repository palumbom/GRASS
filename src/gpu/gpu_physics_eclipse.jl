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

function calc_mu_eclipse_gpu(x, y, z, OPx, OPy, OPz)
    # (x, y, z) is the Sun-to-patch vector and OP the observer-to-patch vector;
    # the cosine of the angle to the observer is minus the cosine between them
    return -calc_mu_gpu(x, y, z, OPx, OPy, OPz)
end

function sky_frame_gpu(OS_bary, sun_rot_mat)
    # unit vectors of the observer's sky frame in the barycentric frame: projected
    # solar north, and west (to the right with north up; the receding limb).
    # OS_bary is observer -> Sun; the third column of the IAU_SUN -> J2000
    # rotation is the solar north pole.
    ux = OS_bary[1]
    uy = OS_bary[2]
    uz = OS_bary[3]
    un = CUDA.sqrt(ux^2.0 + uy^2.0 + uz^2.0)
    ux /= un
    uy /= un
    uz /= un

    px = sun_rot_mat[7]
    py = sun_rot_mat[8]
    pz = sun_rot_mat[9]
    pu = px * ux + py * uy + pz * uz
    nx = px - pu * ux
    ny = py - pu * uy
    nz = pz - pu * uz
    nn = CUDA.sqrt(nx^2.0 + ny^2.0 + nz^2.0)
    nx /= nn
    ny /= nn
    nz /= nn

    # west = line of sight x north
    wx = uy * nz - uz * ny
    wy = uz * nx - ux * nz
    wz = ux * ny - uy * nx
    return nx, ny, nz, wx, wy, wz
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

function sphere_to_cart_gpu_eclipse(ρ, ϕ, θ) #latitude (defined from xy plane rather than from x-axis), longitude
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

function legendreP(l::Int, x::T) where T<:AF
    if l == 0
        return one(T)
    elseif l == 1
        return x
    else
        Pnm2 = one(T)
        Pnm1 = x
        Pn = zero(T)
        for n in 2:l
            Pn = ((2n - 1) * x * Pnm1 - (n - 1) * Pnm2) / n
            Pnm2, Pnm1 = Pnm1, Pn
        end
        return Pn
    end
end

function legendre_dtheta_gpu(l::Int, θ::T) where T<:AF
    # d/dθ of P_l(cos θ) with the normalization of the HMI bulk-velocity fit
    # (Kashyap et al. 2021, arXiv:2105.12055): sqrt(2l+1) / (sqrt(2) sqrt(l(l+1))).
    # θ is the colatitude in radians, measured from the north pole.
    l == 0 && return zero(T)
    cosθ = cos(θ)
    sinθ = sin(θ)

    # dP_l/dz from (1 - z^2) P_l'(z) = l (P_{l-1}(z) - z P_l(z))
    P = legendreP(l, cosθ)
    Plm1 = legendreP(l - 1, cosθ)
    dP = l * (Plm1 - cosθ * P) / (one(T) - cosθ^2)

    # chain rule dz/dθ = -sin θ, then normalize
    norm = sqrt(T(2l + 1)) / (sqrt(T(2)) * sqrt(T(l * (l + 1))))
    return -sinθ * dP * norm
end

function colat_tangent_gpu_eclipse(ϕ, θ) #latitude, longitude; same convention as sphere_to_cart_gpu_eclipse
    # unit vector along increasing colatitude (southward) at the surface point
    sinϕ = sin(ϕ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    cosθ = cos(θ)

    x = sinϕ * cosθ
    y = sinϕ * sinθ
    z = -cosϕ
    return x, y, z
end

function meridional_velocity_gpu(θ, lt, MF1, MF2)
    # line-of-sight meridional flow velocity at colatitude θ (radians). lt is the
    # projection of the southward unit tangent onto the line of sight, positive
    # away from the observer. MF1, MF2 are the s = 2, 4 coefficients of the HMI
    # bulk-velocity fit, in m/s.
    return (MF1 * legendre_dtheta_gpu(2, θ) + MF2 * legendre_dtheta_gpu(4, θ)) * lt
end
