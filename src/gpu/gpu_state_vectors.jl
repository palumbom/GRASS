function calc_state_vector_gpu!(xyz, xyz_dot, epoch, mass_p, period_p, a_p, e_p, i_p, Ω_p, ω_p, L_p)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x

    # get the reduced mass
    μ = mass_p / (mass_p + 1.0)

    # get mean motion
    if CUDA.isnan(period_p)
        n = 2π * CUDA.sqrt(μ / a_p^3.0)
    else
        n = 2π / period_p
    end

    # get longitude of periapsis
    ω̄ = Ω_p + ω_p

    # get the mean anomaly
    M0 = L_p - ω̄

    # loop over time
    for t in idx:sdx:CUDA.length(epoch)
        # get the mean anomaly at epoch
        M = CUDA.mod(M0 + n * epoch[t], 2π)

        # solve kepler's equation for eccentric anomaly
        E = calc_ecc_anom_iterative_laguerre(M, e_p)

        # get the true anomaly
        ν = 2.0 * CUDA.atan(CUDA.sqrt(1.0 + e_p) * CUDA.sin(E / 2.0),
                            CUDA.sqrt(1.0 - e_p) * CUDA.cos(E / 2.0))

        # get distance to central body
        rt = a_p * (1.0 - e_p * CUDA.cos(E))
        rt *= 1.496e11/ 6.957e8

        # assign position in body frame
        @inbounds xyz[1, t] = rt .* CUDA.cos(ν)
        @inbounds xyz[2, t] = rt .* CUDA.sin(ν)
        @inbounds xyz[3, t] = 0.0

        # get velocity magnitude
        # TODO I have no idea if this is right or what units these are
        vt = sqrt(μ * a_p) / rt

        # assign velocity vector
        xyz_dot[1, t] = vt * (-CUDA.sin(E))
        xyz_dot[2, t] = vt * (CUDA.sqrt(1.0 - e_p^2.0) * CUDA.cos(E))
        xyz_dot[3, t] = 0.0

        # shortcut trig evals
        cosω = CUDA.cos(ω_p)
        sinω = CUDA.sin(ω_p)
        cosi = CUDA.cos(i_p)
        sini = CUDA.sin(i_p)
        cosΩ = CUDA.cos(Ω_p)
        sinΩ = CUDA.sin(Ω_p)

        # copy values before rotation
        x0 = CUDA.copy(xyz[1, t])
        y0 = CUDA.copy(xyz[2, t])
        z0 = CUDA.copy(xyz[3, t])

        # rotate into observer frame
        @inbounds xyz[1, t] = x0 * (cosω * cosΩ - sinω * cosi * sinΩ) - y0 * (sinω * cosΩ + cosω * cosi * sinΩ)
        @inbounds xyz[2, t] = x0 * (cosω * sinΩ + sinω * cosi * cosΩ) + y0 * (cosω * cosi * cosΩ - sinω * sinΩ)
        @inbounds xyz[3, t] = x0 * (sinω * sini) + y0 * (cosω * sini)

        # copy values before rotation
        x0 = CUDA.copy(xyz_dot[1, t])
        y0 = CUDA.copy(xyz_dot[2, t])
        z0 = CUDA.copy(xyz_dot[3, t])

        # rotate into observer frame
        @inbounds xyz_dot[1, t] = x0 * (cosω * cosΩ - sinω * cosi * sinΩ) - y0 * (sinω * cosΩ + cosω * cosi * sinΩ)
        @inbounds xyz_dot[2, t] = x0 * (cosω * sinΩ + sinω * cosi * cosΩ) + y0 * (cosω * cosi * cosΩ - sinω * sinΩ)
        @inbounds xyz_dot[3, t] = x0 * (sinω * sini) + y0 * (cosω * sini)
    end
    return nothing
end


function calc_state_vector!(ros_allocs::RossiterAllocsGPU{T1}, body::Planet{T2}) where {T1<:AF, T2<:AF}
    # alias rossiter allocs
    xyz_planet = ros_allocs.xyz_planet
    xyz_dot_planet = ros_allocs.xyz_dot_planet
    xyz_star = ros_allocs.xyz_star
    xyz_dot_star = ros_allocs.xyz_dot_star
    epochs = ros_allocs.epochs

    # get epochs
    @cusync epochs .= CuArray{T1}(collect(range(0.0, length=length(epochs), step = 15.0 / 3.154e7)))

    # get threads and blocks
    threads1 = 256
    blocks1 = cld(length(epochs), prod(threads1))

    # loop over epoch
    @cusync @cuda threads=threads1 blocks=blocks1 calc_state_vector_gpu!(xyz_planet, xyz_dot_planet, epochs,
                                                                         mass(body), period(body), a(body),
                                                                         e(body), i(body), Ω(body),
                                                                         ω(body), L(body))
    return nothing
end
