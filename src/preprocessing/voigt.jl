import SpecialFunctions.erfcx

function voigt(x, sigma, gamma)
    z = (x + im * gamma) / (sqrt(2) * sigma)
    return real(faddeeva(z)) / (sigma * sqrt(2π))
end

function faddeeva(x::Complex{T}) where T<:AF
    return erfcx(-im * x)
end

@. fit_voigt(x, p) = p[1] * voigt(x - p[2], p[3], p[4]) + p[5]
