# set discrete values of mu for input observations
const disc_ax = [:n, :e, :s, :w, :c]

make_grid(N::Integer) = range(-1.0, 1.0, length=N)
make_grid(;N::Integer=256) = make_grid(N)

get_grid_xs(grid::ProductIterator) = getindex.(collect(grid), 1)
get_grid_ys(grid::ProductIterator) = getindex.(collect(grid), 2)

function get_grid_edges(grid::StepRangeLen)
    start = first(grid) - 0.5 * step(grid)
    stop = last(grid) + 0.5 * step(grid)
    return range(start, stop, step=step(grid))
end

function calc_dist2(x1::T, y1::T, x2::T, y2::T) where T<:AF
    return (x1 - x2)^2 + (y1 - y2)^2
end

function calc_dist2(t1::Tuple{T,T}, t2::Tuple{T,T}) where T<:AF
    return calc_dist2(t1..., t2...)
end

function calc_r2(x::T,y::T) where T<:AF
    return x*x + y*y
end

function calc_r2(t::Tuple{T,T}) where T<:AF
    return calc_r2(t[1], t[2])
end

# Get mu on disk from provided (x,y) position:
# TODO: relies on small angle approx?
function calc_mu(r2::T) where T<:AF
    return sqrt((r2 < one(T)) * (one(T) - r2))
end

function calc_mu(x::T, y::T) where T<:AF
    return calc_mu(calc_r2(x,y))
end

function calc_mu(t::Tuple{T,T}) where T<:AF
    return calc_mu(calc_r2(t))
end

# Calculate mu for each position on a grid
function mu_map(grid::AA{T,1}) where T<:AF
    return calc_mu.((x,y) for x in grid, y in grid)
end

function mu_map(N::Integer)
    return mu_map(make_grid(N))
end

function calc_norm_term(μ::T, N::Integer, u1::T, u2::T) where T<:AF
    return quad_limb_darkening(μ,u1,u2) * π / (2.0 * N^2)
end

function calc_norm_term(x::T, y::T, N::Integer, u1::T, u2::T) where T<:AF
    return calc_norm_term(calc_mu(x,y), N, u1, u2)
end

function calc_norm_term(x::T, y::T, disk::DiskParams) where T<:AF
    return calc_norm_term(calc_mu(x,y), disk.N, disk.u1, disk.u2)
end


function calc_norm_terms(N::Int, u1::T, u2::T; nsubgrid::Int=256) where T<:AF
    # get grids
    len = Int(N/2)
    grid = Iterators.product(range(0.0, 1.0, length=len),
                             range(0.0, 1.0, length=len))
    grid_range = range(0.0, 1.0, length=len)
    grid_edges = get_grid_edges(grid_range)

    # allocate memory for norm_terms
    norm_terms = zeros(N, N)
    norm_terms_quadrant = zeros(size(grid))

    # calculate distances
    dist2 = map(z -> calc_dist2(z, (0.0, 0.0)), grid)

    # find edge pixels
    is_edge = (dist2 .>= (1.0 - step(grid_range))) .&& (dist2 .<= 1.0 + step(grid_range))
    is_not_edge = (.!is_edge) .&& (dist2 .< 1.0)

    # calculate norm terms on not_edge pixels
    for idx in findall(is_not_edge) #CartesianIndices(norm_terms)
        i = idx[1]; j= idx[2]
        norm_terms_quadrant[idx] = quad_limb_darkening.(grid_range[i], grid_range[j], u1, u2)
    end

    # calculate norm terms for edge pixels
    for idx in findall(is_edge)
        i = idx[1]; j= idx[2]

        # make subgrid
        xrange = range(grid_edges[i], grid_edges[i+1], length=nsubgrid)
        yrange = range(grid_edges[j], grid_edges[j+1], length=nsubgrid)
        subgrid = Iterators.product(xrange, yrange)

        # calculate norm terms on subgrid
        sublimbdarks = quad_limb_darkening.(calc_mu.(subgrid), u1, u2)
        norm_terms_quadrant[idx] = mean(sublimbdarks)
    end

    # now copy quadrants over
    norm_terms[len+1:end, len+1:end] .= norm_terms_quadrant
    norm_terms[1:len, len+1:end] .= reverse!(norm_terms_quadrant, dims=(1))
    norm_terms[len+1:end, 1:len] .= reverse!(norm_terms_quadrant, dims=(1,2))
    norm_terms[1:len, 1:len] .= reverse!(norm_terms_quadrant, dims=(1))
    return norm_terms ./ sum(norm_terms)
end

function calc_norm_terms(disk::DiskParams; nsubgrid::Int=256)
    return calc_norm_terms(disk.N, disk.u1, disk.u2, nsubgrid=nsubgrid)
end


# Find the nearest axis to a given point on a grid
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

function find_nearest_mu(mu::T, disc_mu::AA{T,1}) where T<:AF
    return searchsortednearest(disc_mu, mu)
end

function find_nearest_mu(x::T, y::T, disc_mu::AA{T,1}) where T<:AF
    return find_nearest_mu(calc_mu(x,y), disc_mu)
end

function assemble_dict_key(mu_ind::Int, ax::Symbol, mu_symb::AA{Symbol,1})
    if mu_symb[mu_ind] == :mu10 # e.g., if mu == 1.0
        return (:c, :mu10)
    else
        return (ax, mu_symb[mu_ind])
    end
end

function get_key_for_pos(x::T, y::T, disc_mu::AA{T,1}, mu_symb::AA{Symbol,1}) where T<:AF
    # make sure we are not off the disk
    if x^2 + y^2 > 1
        return nothing
    end

    # find nearest mu index
    mu_ind = find_nearest_mu(x, y, disc_mu)

    # get nearest axis and mu
    near_mu = disc_mu[mu_ind]
    near_ax = find_nearest_ax(x, y)

    # assemble dictionary key & extract values
    key = assemble_dict_key(mu_ind, near_ax, mu_symb)
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
