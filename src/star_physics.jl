function absorption_line(x::T; continuum=one(T), mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    return continuum - depth * exp(-(x-mid)^2.0/(2.0 * width^2.0))
end

function gaussian_line(x::T; continuum=one(T), mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    return absorption_line(x, continuum=continuum, mid=mid, width=width, depth=depth)
end

function lorentzian_line(x::T; mid=zero(T), width=one(T), depth=one(T)) where T<:AF
    y = (x - mid)/(width/2)
    return one(T) - (depth / (one(T) + y^2))
end

# calculate doppler factor
function doppler_factor(vel::T) where T<:AF
    num = one(T) + vel/c_ms
    den = one(T) - vel/c_ms
    return sqrt(num/den)
end

# Width of line due to thermal broadening
function width_thermal(; λ::T1=1.0, M::T1=1.0, T::T1=5778.0, v_turb::T1=0.0) where T1<:AF
    return sqrt((2.0*kB/mH) * (T/M) + v_turb^2) * (λ/c)
end

function quad_limb_darkening_eclipse(μ::T, wavelength::T) where T<:AF
    μ < zero(T) && return 0.0    

    # index = findmin(x->abs(x-wavelength), lambda_nm)[2]

    # return a0[index] + a1[index]*μ + a2[index]*μ^2 + a3[index]*μ^3 + a4[index]*μ^4 + a5[index]*μ^5
    
    index = findmin(x->abs(x-wavelength), Kostogryz_LD_file[!, "wavelength"])[2]
    Kostogryz_LD_array = Kostogryz_LD_file[index, ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]]
    Kostogryz_LD_interpol = linear_interp([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], collect(Kostogryz_LD_array))

    return Kostogryz_LD_interpol(μ)
end

function quad_limb_darkening_NIR(μ::T) where T
    """
    limb darkening prescription for NIR based on mu angle  
    """
    μ < zero(T) && return 0.0
    return 0.59045 + 1.41938*μ - 3.01866*μ^2 + 3.99843*μ^3 - 2.67727*μ^4 + 0.068758*μ^5
end

function quad_limb_darkening(μ::T, u1::T, u2::T) where T<:AF
    μ < zero(T) && return 0.0
    return !iszero(μ) * (one(T) - u1*(one(T)-μ) - u2*(one(T)-μ)^2)
end

"""
    rotation_period(sin_lat; A=14.713, B=2.396, C=1.787)

Calculate stellar rotation period in days at given sine of latitude.

# Arguments
- `sin_lat::Float64`: Sine of stellar latitude.
- `A::Float64=14.713`: Coefficient from Snodgrass & Ulrich (1990)
- `B::Float64=2.396`: Coefficient from Snodgrass & Ulrich (1990)
- `C::Float64=1.787`: Coefficient from Snodgrass & Ulrich (1990)
"""
function rotation_period(ϕ::T; A::T=14.713, B::T=-2.396, C::T=-1.787) where T<:AF
    @assert -π/2 <= ϕ <= π/2
    sinϕ = sin(ϕ)
    return 360.0/(A + B * sinϕ^2.0 + C * sinϕ^4.0)
end

function v_scalar(lat, lon)
    return (2π * sun_radius * cos(lat)) / rotation_period(lat)
end

function projected!(A::Matrix, B::Matrix, out_no_cb::Matrix)
    """
    determine projected velocity of each cell onto line of sight to observer - serial

    A: matrix with xyz and velocity of each cell
    B: matrix with line of sight from each cell to observer
    out: matrix of projected velocities
    """
    for i in 1:length(A)
        vel = A[i][4:6]
        angle = dot(B[i][1:3], vel) / (norm(B[i][1:3]) * norm(vel))

        out_no_cb[i] = (norm(vel) * angle)
    end
    return
end


"""
    patch_velocity_los(x, y; rstar, pole)

Compute line of sight velocity of a patch of stellar surface given by x,y (assumed in [-1,1]).
Return value is in (Rsol/day)/speed of light (i.e., dimensionless like z = v/c)

# Arguments
- `rstar::Float64=1.0`: in R_sol (affects return velocity, but not x,y)
- `pole = (0,1,0)`: unit vector for stellar rotation axis (default is equator-on)
"""
function patch_velocity_los(ϕ::T, θ::T, disk::DiskParams{T}; P⃗=[0.0, disk.ρs, 0.0]) where T<:AF
    # get vector pointing from spherical circle to patch
    xyz = sphere_to_cart(disk.ρs, ϕ, θ)

    # get vector from spherical circle center to surface patch
    C⃗ = xyz .- [0.0, xyz[2], 0.0]

    # get magnitude of velocity vector
    v0 = 2π * disk.ρs * cos(ϕ) / rotation_period(ϕ; A=disk.A, B=disk.B, C=disk.C)

    # get in units of c
    v0 /= c_Rsun_day

    # get velocity vector direction and set magnitude
    vel = cross(C⃗, P⃗)
    vel /= norm(vel)
    vel *= v0

    # rotate by stellar inclination
    xyz .= disk.R_x * xyz
    vel .= disk.R_x * vel

    # find get vector from observer to surface patch, return projection
    O⃗_surf = xyz .- disk.O⃗
    angle = dot(O⃗_surf, vel) / (norm(O⃗_surf) * norm(vel))
    return norm(vel) * angle
end