function calc_mu_gpu(xyz, O⃗)
    dp = xyz[1] * O⃗[1] + xyz[2] * O⃗[2] + xyz[3] * O⃗[3]
    n1 = CUDA.sqrt(O⃗[1]^2.0 + O⃗[2]^2.0 + O⃗[3]^2.0)
    n2 = CUDA.sqrt(xyz[1]^2.0 + xyz[2]^2.0 + xyz[3]^2.0)
    return dp / (n1 * n2)
end

function sphere_to_cart_gpu!(xyz, ρs, ϕc, θc)
    # compute trig quantities
    sinϕ = CUDA.sin(ϕc)
    sinθ = CUDA.sin(θc)
    cosϕ = CUDA.cos(ϕc)
    cosθ = CUDA.cos(θc)

    # now get cartesian coords
    @inbounds xyz[1] = ρs * cosϕ * cosθ
    @inbounds xyz[2] = ρs * cosϕ * sinθ
    @inbounds xyz[3] = ρs * sinϕ
    return nothing
end

function rotate_vector_gpu!(xyz, R_θ)
    # parse out components
    x = xyz[1]
    y = xyz[2]
    z = xyz[3]

    # do dot product
    @inbounds xyz[1] = x * R_θ[1,1] + y * R_θ[1,2] + z * R_θ[1,3]
    @inbounds xyz[2] = x * R_θ[2,1] + y * R_θ[2,2] + z * R_θ[2,3]
    @inbounds xyz[3] = x * R_θ[3,1] + y * R_θ[3,2] + z * R_θ[3,3]
    return nothing
end


function rotation_period_gpu(ϕ, A, B, C)
    sinϕ = sin(ϕ)
    return 360.0/(A + B * sinϕ^2.0 + C * sinϕ^4.0)
end
