function calc_mu_gpu(x, y, z, O⃗)
    dp = x * O⃗[1] + y * O⃗[2] + z * O⃗[3]
    n1 = CUDA.sqrt(O⃗[1]^2.0 + O⃗[2]^2.0 + O⃗[3]^2.0)
    n2 = CUDA.sqrt(x^2.0 + y^2.0 + z^2.0)
    return dp / (n1 * n2)
end

function sphere_to_cart_gpu(ρs, ϕ, θ)
    # compute trig quantities
    sinϕ = CUDA.sin(ϕ)
    sinθ = CUDA.sin(θ)
    cosϕ = CUDA.cos(ϕ)
    cosθ = CUDA.cos(θ)

    # now get cartesian coords
    x = ρs * cosϕ * cosθ
    y = ρs * cosϕ * sinθ
    z = ρs * sinϕ
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
