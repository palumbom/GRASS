# set discrete values of mu for input observations
const disc_ax = [:n, :e, :s, :w, :c]

# make_grid(N::Integer) = range(-1.0, 1.0, length=N)
# make_grid(;N::Integer=256) = make_grid(N)

function make_grid(N::Integer)
    # create grid edges
    ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=N)
    θe = range(deg2rad(0.0), deg2rad(360.0), length=N)
    return ϕe, θe
end

function get_grid_centers(grid::StepRangeLen)
    start = first(grid) + 0.5 * step(grid)
    stop = last(grid) - 0.5 * step(grid)
    return range(start, stop, step=step(grid))
end

function calc_area_element(ρs::T, ϕc::T, dϕ::T, dθ::T) where T<:AF
    return ρs^2.0 * sin(π/2.0 - ϕc) * dϕ * dθ
end

function calc_projected_area_element(ϕc::T, θc::T, disk::DiskParams{T}) where T<:AF
    # get area element
    dA = calc_area_element(disk.ρs, ϕc, step(disk.ϕe), step(disk.θe))

    # get cartesian coords and rotate them
    xyz = sphere_to_cart(disk.ρs, ϕc, θc)
    xyz .= disk.R_θ * xyz

    # get vector from observer to surface element and return projection
    O⃗_surf = xyz .- disk.O⃗
    return dA * abs(dot(O⃗_surf, xyz))
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

function calc_mu(ϕ::T, θ::T; R_θ::AA{T,2}=Matrix(1.0I,3,3), O⃗::AA{T,1}=[0.0, 220.0, 0.0]) where T<:AF
    # get cartesian coords and rotate them
    xyz = sphere_to_cart(one(T), ϕ, θ)
    xyz .= R_θ * xyz
    return dot(O⃗, xyz) / (norm(O⃗) * norm(xyz))
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

function find_nearest_ax_code(x::T, y::T) where T<:AF
    if (x^2 + y^2) > one(T)
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

function get_key_for_pos(μ::T, ϕ::T, θ::T, disc_mu::AA{T,1}, disc_ax::AA{Int,1}; R_θ::AA{T,2}=Matrix(1.0I,3,3),) where T<:AF
    # make sure we are not off the disk
    if μ < 0.0
        return nothing
    end

    # get cartesian position for axis
    xyz = sphere_to_cart(one(T), ϕ, θ)
    xyz .= R_θ * xyz

    # find the nearest mu ind and ax code
    mu_ind = searchsortednearest(disc_mu, μ)
    ax_val = find_nearest_ax_code(xyz[1], xyz[2])

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
