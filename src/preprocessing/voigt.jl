import SpecialFunctions.erfcx

function gaussian(x, amp, mu, sigma)
    return amp * exp(-0.5 * ((x - mu)/sigma)^2)
end


function lorentzian(x, amp, mu, gamma)
    return amp * (gamma^2 / ((x - mu)^2 + gamma^2))
end

function faddeeva(x::Complex{T}) where T<:AF
    return erfcx(-im * x)
end

function voigt(x, amp, mu, sigma, gamma)
    if iszero(sigma)
        return lorentzian(x, amp, mu, gamma)
    elseif iszero(gamma)
        return gaussian(x, amp, mu, sigma)
    else
        z = ((x - mu) + im * gamma) / (sqrt(2) * sigma)
        return amp * real(faddeeva(z)) / (sigma * sqrt(2π))
    end
    return nothing
end

@. fit_voigt(x, p) = exp(-voigt(x, p[1], p[2], p[3], p[4]))
