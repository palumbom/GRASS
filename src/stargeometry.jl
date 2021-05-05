"""
Author: Michael Palumbo
Created: May 2019
Contact: mlp95@psu.edu
"""

# set discrete values of mu for input observations
const disc_mu = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0]
const mu_symb = [:mu02, :mu03, :mu04, :mu05, :mu06, :mu07, :mu08, :mu085, :mu09, :mu095, :mu10]
const disc_ax = [:n, :e, :s, :w, :c]

"""
Make iterator representing uniform NxN grid
"""
make_grid(N::Integer) = range(-1.0, 1.0, length=N)
make_grid(;N::Integer=256) = make_grid(N)

"""
Calculate length of distance vector projected on face of star
"""
function calc_r2(x::T,y::T) where T<:AF
    return x*x + y*y
end

function calc_r2(t::Tuple{T,T}) where T<:AF
    return calc_r2(t[1], t[2])
end


"""
Get mu on disk from provided (x,y) position:
(x,y) = position on projected disk
TODO: relies on small angle approx
"""
function calc_mu(t::Tuple{T,T}) where T<:AF
    r2 = calc_r2(t)
    return sqrt((r2 < one(T)) * (one(T) - r2))
end

function calc_mu(x::T, y::T) where T<:AF
    return calc_mu((x,y))
end

"""
Calculate mu for each position on a grid
"""
function mu_map(grid::AA{T,1}) where T<:AF
    return calc_mu.((x,y) for x in grid, y in grid)
end

function mu_map(N::Integer)
    return mu_map(make_grid(N))
end

"""

"""
function norm_term(x::T, y::T, N::Integer, u1::T, u2::T) where T<:AF
    return quad_limb_darkening(x,y,u1,u2) * π / (2.0 * N^2)
end

function norm_term(x::T, y::T, disk::DiskParams) where T<:AF
    return norm_term(x, y, disk.N, disk.u1, disk.u2)
end

"""
Find the nearest axis to a given point on a grid
"""
function find_nearest_ax(x::T, y::T) where T<:AF
    if (x^2 + y^2) > one(T)
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

function find_nearest_mu(x::T, y::T) where T<:AF
    return searchsortednearest(disc_mu, calc_mu(x, y))
end

function assemble_dict_key(mu_ind::Int, ax::Symbol; mu_symb::AA{Symbol,1}=mu_symb)
    if mu_ind == 11 # e.g., if mu == 1.0
        return (:c, :mu10)
    else
        return (ax, mu_symb[mu_ind])
    end
end

function get_key_for_pos(x::T, y::T) where T<:AF
    # find nearest mu index
    mu_ind = find_nearest_mu(x, y)

    # get nearest axis and mu
    near_mu = disc_mu[mu_ind]
    near_ax = find_nearest_ax(x, y)

    # assemble dictionary key & extract values
    key = assemble_dict_key(mu_ind, near_ax)
    return key
end

function mu_to_xy(mu::T, ax::Symbol) where T<:AF
    @assert zero(T) <= mu <= one(T)
    @assert ax in [:n, :e, :s, :w, :c]

    # get radial distance from disk center
    r = sqrt(one(T) - mu^2)


    if ax == :n
        return (0.0, r)
    elseif ax == :e
        return (-r, 0.0)
    elseif ax == :w
        return (r, 0.0)
    elseif ax == :s
        return (0.0, -r)
    else
        return (0.0, 0.0)
    end
end
