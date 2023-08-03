# set discrete values of mu for input observations
const disc_ax = [:n, :e, :s, :w, :c]

function make_grid(N::Integer)
    # create grid edges
    ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=N)
    θe = range(deg2rad(0.0), deg2rad(360.0), length=N)
    return ϕe, θe
end

function get_grid_centers(grid::StepRangeLen)
    start = first(grid) + 0.5 * step(grid)
    stop = last(grid) - 0.5 * step(grid)
    return range(start, stop, length=length(grid)-1)
end

function calc_dA(ρs::T, ϕc::T, dϕ::T, dθ::T) where T<:AF
    return ρs^2.0 * sin(π/2.0 - ϕc) * dϕ * dθ
end

function calc_projected_dA(ϕc::T, θc::T, disk::DiskParams{T}) where T<:AF
    # get area element
    dA = calc_dA(disk.ρs, ϕc, step(disk.ϕe), step(disk.θe))

    # get cartesian coords and rotate them
    xyz = sphere_to_cart(disk.ρs, ϕc, θc)
    xyz .= disk.R_θ * xyz

    # get vector from observer to surface element and return projection
    O⃗_surf = xyz .- disk.O⃗
    return dA * abs(dot(O⃗_surf, xyz))
end

function sphere_to_cart(ρ::T, ϕ::T, θ::T) where T
    # compute trig quantitites
    sinϕ = sin(ϕ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    cosθ = cos(θ)

    # now get cartesian coords
    x = ρ * cosϕ * cosθ
    y = ρ * cosϕ * sinθ
    z = ρ * sinϕ
    return [x, y, z]
end

function calc_mu(xyz::AA{T,1}, R_θ::AA{T,2}, O⃗::AA{T,1}) where T<:AF
    # get cartesian coords and rotate them
    xyz .= R_θ * xyz
    return dot(O⃗, xyz) / (norm(O⃗) * norm(xyz))
end

# Find the nearest axis to a given point on a grid
function find_nearest_ax(x::T, y::T) where T<:AF
    if (x^2.0 + y^2.0) > one(T)
        return :off
    elseif ((y == zero(T)) & (x == zero(T))) # center
        return :c
    elseif y >= abs(x) # north
        return :n
    elseif y <= -abs(x) # south
        return :s
    elseif x <= -abs(y) # east
        return :e
    elseif x >= abs(y) # west
        return :w
    end
end

function find_nearest_ax_code(x::T, y::T) where T<:AF
    if (x^2.0 + y^2.0) > one(T)
        return nothing
    elseif ((y == zero(T)) & (x == zero(T))) # center
        return 0
    elseif y >= abs(x) # north
        return 1
    elseif y <= -abs(x) # south
        return 2
    elseif x <= -abs(y) # east
        return 3
    elseif x >= abs(y) # west
        return 4
    end
end

function ax_code_to_symbol(code::Int)
    return [:c, :n, :s, :e, :w][code+1]
end

function get_key_for_pos(μ::T, x::T, y::T, disc_mu::AA{T,1}, disc_ax::AA{Int,1}) where T<:AF
    # make sure we are not off the disk
    if μ <= 0.0
        return (:off, :off)
    end

    # find the nearest mu ind and ax code
    mu_ind = searchsortednearest(disc_mu, μ)
    ax_val = find_nearest_ax_code(x, y)

    # return early if the nearest mu is 1.0
    if disc_mu[mu_ind] == 1.0
        return (:c, :mu10)
    end

    # find subarray of disc_mu and disk_ax matching mu
    idxs = findall(disc_mu .== disc_mu[mu_ind])
    mu_view = view(disc_mu, idxs)
    ax_view = view(disc_ax, idxs)

    # move to new axis if it isn't present in the data
    if !(ax_val in ax_view)
        ax_val = ax_view[1]
    end

    # convert mu and ax codes to symbol key
    mu_symb = mu_to_symb(disc_mu[mu_ind])
    ax_symb = ax_code_to_symbol(ax_val)
    return (ax_symb, mu_symb)
end

function key_to_code(key, disc_mu, disc_ax)
    ax_sym = key[1]
    mu_sym = key[2]

    ax_code = parse_ax_string(ax_sym)
    mu_code = parse_mu_string(mu_sym)

    idxs1 = findall(ax_code .== disc_ax)
    idxs2 = findall(mu_code .== disc_mu)
    return collect(intersect(Set(idxs1), Set(idxs2)))[1]
end
