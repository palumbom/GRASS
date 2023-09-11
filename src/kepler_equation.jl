"""
Function for numerically solving Kepler's Equation: M = E - e * sin(E)

Meant to be passed to rootsolver function

# Arguments
- `E::Float64`: the eccentric anomaly in radians (to be found via root solving).
- `e::Float64`: the eccentricity.
- `M::Float64`: the mean anomaly in radians.
"""
function kepler(E::Float64, e::Float64, M::Float64)
    @assert 0.0 <= e < 1.0
    return E - e * sin(E) - M
end

"""
Compute an initial guess for eccentric anomaly.
Called only by `calc_ecc_anom_iterative_laguerre()`.
(Based on "The Solution of Kepler's Equations - Part Three",
Danby, J. M. A. (1987) Journal: Celestial Mechanics,
Volume 40, Issue 3-4, pp. 303-312, 1987CeMec..40..303D)

# Arguments
- `M::Float64`: the mean anomaly in radians
- `e::Float64`: the eccentricity
"""
function ecc_anom_init_guess_danby(M::Float64, e::Float64)
    if M < 0.1
        return M + ((6.0 * M)^(1.0/3.0) - M) * e^2
    else
        return M + 0.85 * e
    end
end

"""
Update the current guess for the solution to Kepler's equation
Called only by `calc_ecc_anom_iterative_laguerre()`.
(Based on "An Improved Algorithm due to Laguerre for the Solution of Kepler's Equation",
Conway, B. A.  (1986) Celestial Mechanics,
Volume 39, Issue 2, pp.199-211, 1986CeMec..39..199C)

# Arguments
- `E::Float64`: the eccentric anomaly in radians
- `M::Float64`: the mean anomaly in radians
- `e::Float64`: the eccentricity
"""
function update_ecc_anom_laguerre(E::Float64, M::Float64, e::Float64)
    esin = e * sin(E)
    ecos = e * cos(E)
    F = (E - esin) - M
    Fp = 1.0 - ecos
    Fpp = esin
    n = 5
    root = sqrt(abs((n-1)*((n-1)*Fp*Fp-n*F*Fpp)))

    if Fp > 0.0
        denom = Fp + root
    else
        denom = Fp - root
    end

    return E - n*F/denom
end

"""
Iteratively calculates a value for eccentric anomaly, converges within a given tolerance.

# Arguments
- `mean_anom::Float64`: the mean anomaly in radians
- `ecc::Float64`: the eccentricity

# Keyword Arguments
- `tol::Float64`: the tolerance for convergence, default 1.0e-8

# Output
- `E::Float64`: the eccentric anomaly in radians
"""
function calc_ecc_anom_iterative_laguerre(mean_anom::Float64, ecc::Float64, tol::Float64=1.0e-8, maxit::Int64=200)
    M = mod2pi(mean_anom)
    E = ecc_anom_init_guess_danby(M, ecc)
    E_old = E
    for i in 1:maxit
        E_old = E
        E = update_ecc_anom_laguerre(E_old, M, ecc)

        # end loop if tolerance achieved
        if abs(E - E_old) < tol break end
    end

    # toss error if convergence not achieved
    @assert abs(E - E_old) < tol
    return E
end
