function step_forward_matrix(omega::T, gamma::T, dt::T) where T<:AF
    # return [x, rv] after time dt
    Minv = inv([1.0 0.0; -0.5 * gamma omega])

    A = zeros(2,2)
    A[1] = cos(omega * dt)
    A[2] = -omega * sin(omega * dt) - 0.5 * gamma * cos(omega * dt)
    A[3] = sin(omega * dt)
    A[4] = omega * cos(omega * dt) - 0.5 * gamma * sin(omega * dt)

    A .*= exp(-0.5 * gamma * dt)
    return A * Minv
end

function take_one_step!(x_v, step_half_matrix, step_full_matrix, driving_amp)
    # returns x0.5, v0.5 (prediction for next observation time)
    # and x1, v1 (starting position post-kick)
    x_half, v_half = step_half_matrix * x_v
    x_one, v_one = step_full_matrix * x_v
    v_one += driving_amp * randn()

    # set vector elements
    x_v[1] = x_one
    x_v[2] = v_one
    return x_half, v_half
end

function take_many_steps(omega, gamma, dt, driving_amp; timescale=365.0)
    # timescale is number of days to observe
    # assert dt < np.pi / omega, "ERROR: you're not well-sampled. decrease dt."
    # if dt > 0.5 * np.pi / omega: print("WARNING: your coarse time spacing makes even cubic spline risky")

    ts = (dt .* range(1.0, timescale * 24. * 3600. / dt)) .- dt
    xs = zeros(eltype(ts), size(ts)...)
    rvs = zeros(eltype(ts), size(ts)...)
    
    step_half_matrix = step_forward_matrix(omega, gamma, 0.5 * dt)
    step_full_matrix = step_forward_matrix(omega, gamma, dt)

    x_v = zeros(2)
    for i in eachindex(ts)
        xs[i], rvs[i] = take_one_step!(x_v, step_half_matrix, step_full_matrix, driving_amp)
    end
    return ts, xs, rvs
end

function simulate_exposure(ts, rvs, start_time, exp_time)
    pad = 100. # seconds - ARBITRARY
    smaller_inds = (ts .> (start_time - pad)) .& (ts .< (start_time + exp_time + pad))
    itp = Spline1D(ts[smaller_inds], rvs[smaller_inds]; w=ones(length(ts[smaller_inds])), k=3, bc="nearest")    
    tiny = 0.1 # 100 ms
    fine_ts = range(start_time, start_time+exp_time, step=tiny) # fine grid
    fine_rvs = itp.(fine_ts)
    return sum(fine_rvs) / length(fine_rvs) # ASSUMES EVEN WEIGHTING - technically incorrect for last point
end