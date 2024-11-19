struct PModes{T<:AF}
    a_max::T
    ν_max::T
    c_env::T
    Δν::T
    γ::T
    f_grid::AA{T,1}
    ω_grid::AA{T,1}
    a_grid::AA{T,1}
    driving_amp_grid::AA{T,1}
    ts::AA{T,1}
    rvs::AA{T,1}
end

function PModes(;a_max=0.55, ν_max=3.1e-3, c_env=0.331e-3, Δν=0.00013, γ=0.0)
    # hard code dt and γ for now
    dt = 15.0
    γ = 1.0 / (2.0 * 24. * 60.0 * 60.0)

    # grid of frequencies, etc.
    f_grid = range(ν_max - 0.001, ν_max + 0.001, step=Δν) 
    l_grid = vcat(zeros(size(f_grid)), ones(size(f_grid)))
    f_grid = vcat(f_grid, f_grid)
    f_grid .+= l_grid .* 0.5 .* Δν # Hogg!
    ω_grid = 2.0 .* π .* f_grid # angular frequencies
    a_grid = a_max^2.0 .* exp.(.-(f_grid .- ν_max) .^ 2.0 ./ (2.0 .* c_env .^2.0)) # amplitudes in m/s
    a_grid .-= 0.4 .* a_grid .* l_grid # Hogg!
    driving_amp_grid = sqrt.(a_grid .* gamma .* dt)

    # hardcode a timespan (in days)
    timespan = 10.0
    ts, xs, rvs = take_many_steps(ω_grid, gamma, dt, driving_amp_grid, timescale=timespan)

    return PModes(a_max, ν_max, c_env, Δν, γ, f_grid, ω_grid, a_grid, driving_amp_grid, ts, rvs)
end
