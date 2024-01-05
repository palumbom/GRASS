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

function calc_proj_dist2(p1::AA{T,1}, p2::AA{T,1}) where T<:AF
    x1 = p1[1]
    x2 = p2[1]
    y1 = p1[2]
    y2 = p2[2]
    return (x1 - x2)^2.0 + (y1-y2)^2.0
end

function calc_dA(ρs::T, ϕc::T, dϕ::T, dθ::T) where T<:AF
    return ρs^2.0 * sin(π/2.0 - ϕc) * dϕ * dθ
end

function sphere_to_cart_eclipse(ρ::T, ϕ::T, θ::T) where T
    # compute trig quantitites
    sinϕ, cosϕ = sincos(ϕ)
    sinθ, cosθ = sincos(θ)

    # now get cartesian coords
    y = ρ * cosϕ * sinθ
    z = ρ * sinϕ
    x = ρ * cosϕ * cosθ
    return [x, y, z]
end

function pole_vector_grid!(A::Matrix, out::Matrix)
    """
    remove the z component of each cell - serial

    A: matrix of xyz orientation of each cell
    out: matrix of xyz orientation with z removed
    """ 
    for i in 1:length(A)
        out[i] = A[i] - [0.0, 0.0, A[i][3]]
    end
    return
end  

function v_vector(A::Matrix, B::Matrix, C::Matrix, out::Matrix)
    """
    determine velocity vector (direction and magnitude) of each cell - serial 

    A: xyz position of cell
    B: xyz position of cell with z removed
    C: scalar velocity of each cell
    out: matrix with xyz and velocity of each cell
    """
    for i in eachindex(A)
        cross_product = cross([0.0,0.0,sun_radius], B[i])
        cross_product ./= norm(cross_product)
        cross_product .*= C[i]
        out[i] = cross_product
    end
    return
end

function sphere_to_cart(ρ::T, ϕ::T, θ::T) where T
    # compute trig quantitites
    sinϕ = sin(ϕ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    cosθ = cos(θ)

    # now get cartesian coords
    x = ρ * cosϕ * sinθ
    y = ρ * sinϕ
    z = ρ * cosϕ * cosθ
    return [x, y, z]
end

function calc_mu(xyz::AA{T,1}, O⃗::AA{T,1}) where T
    return dot(O⃗, xyz) / (norm(O⃗) * norm(xyz))
end

function calc_mu_grid!(A::Matrix, B::Matrix, out::Matrix)
    """
    create matrix of mu values for each cell - serial

    A: matrix of vectors from sun center to cell
    B: matrix of vectors from observer to cell
    out: mu value between A and B
    """
    for i in eachindex(A)
        out[i] = calc_mu(view(A[i], 1:3), view(B[i], 1:3))
        end
    return
end

function find_nearest_ax_code_eclipse(y::T, z::T) where T<:AF
    # if (y^2.0 + z^2.0) > one(T) 
    #     return nothing
    if ((z == zero(T)) & (y == zero(T))) # center
        return 0
    elseif z >= abs(y) # north
        return 1
    elseif z <= -abs(y) # south
        return 2
    elseif y <= -abs(z) # east
        return 3
    elseif y >= abs(z) # west
        return 4
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

function get_key_for_pos(μ::T, ax::Int, disc_mu::AA{T,1}, disc_ax::AA{Int,1}) where T<:AF
    # make sure we are not off the disk
    if μ <= 0.0
        return (:off, :off)
    end

    # find the nearest mu and return early if near disk center
    mu_ind = searchsortednearest(disc_mu, μ)
    if disc_mu[mu_ind] == 1.0
        return (:c, :mu10)
    end

    # copy ax code so we don't mutate it in host array somehow
    ax_val = copy(ax)

    # find subarray of disc_mu and disk_ax matching mu
    idxs = findall(disc_mu .== disc_mu[mu_ind])
    mu_view = view(disc_mu, idxs)
    ax_view = view(disc_ax, idxs)

    # move to new axis if it isn't present in the data
    if !(ax_val in ax_view)
        ax_val = first(ax_view)
    end

    # convert mu and ax codes to symbol key
    mu_symb = mu_to_symb(disc_mu[mu_ind])
    ax_symb = ax_code_to_symbol(ax_val)
    return (ax_symb, mu_symb)
end

function key_to_code(key, soldata)
    disc_mu = soldata.mu
    disc_ax = soldata.ax

    if key == (:off, :off)
        return 0
    end
    ax_sym = key[1]
    mu_sym = key[2]

    ax_code = parse_ax_string(ax_sym)
    mu_code = parse_mu_string(mu_sym)

    idxs1 = findall(ax_code .== disc_ax)
    idxs2 = findall(mu_code .== disc_mu)
    return collect(intersect(Set(idxs1), Set(idxs2)))[1]
end

function calc_proj_dist(p1, p2)
    return acos(dot(p1, p2)/(norm(p1)*norm(p2)))
end