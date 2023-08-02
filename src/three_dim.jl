using LinearAlgebra
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()

function get_xyz_for_surface(ρ::T, ϕ::T, θ::T) where T
    # pre-compute trig quantitites
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

function differential_rotation(ρ::T, ϕ::T) where T
    # get trig quantities
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)

    # get rotation period in days at specified latitude
    A = 14.713
    B = -2.396
    C = -1.787
    period = 360.0/(A + B * sinϕ^2 + C * sinϕ^4)

    # get total velocity vector length
    return 0.000168710673 * 3e8 / period
end

# set the inclination of the star and get rotation matrix
ρₛ = 1.0
iₛ = deg2rad(45.0)
M = [1.0 0.0 0.0;
     0.0 cos(iₛ) sin(iₛ);
     0.0 -sin(iₛ) cos(iₛ)]

# set up lat long grid edges
ϕs = deg2rad.(range(-90.0, 90.0, length=181))
θs = deg2rad.(range(0.0, 360.0, length=361))

# get grid centers
θc = (θs[2:end] + θs[1:end-1])/2.0
ϕc = (ϕs[2:end] + ϕs[1:end-1])/2.0

# get difference elements
dθ = diff(θs)
dϕ = diff(ϕs)

# pole vector in star frame
pole_vec = [0.0, 0.0, ρₛ]

# approx earth-sun distance in solar radii
obs = [0.0, 220.0, 0.0]

# get velocity vector in star frame
μs = zeros(length(ϕc), length(θc))
vs = similar(μs)
dA_proj = similar(vs)

for i in 1:length(ϕc)
    for j in 1:length(θc)
        # get cartesian coords
        xyz = get_xyz_for_surface(ρₛ, ϕc[i], θc[j])

        # get vector pointing from star origin to annulus height
        ann_vec = xyz .- [0.0, 0.0, last(xyz)]

        # take cross product to get vector in direction of surface rotation
        vel = cross(ann_vec, pole_vec)
        vel /= norm(vel)
        vel *= differential_rotation(ρₛ, ϕs[i])

        # rotate by stellar inclination
        xyz_new = M * xyz
        vel_new = M * vel

        # get mu
        μs[i,j] = dot(obs, xyz_new) / (norm(obs) * norm(xyz_new))

        # find get vector from observer to surface patch
        v_obs_surf = xyz_new .- obs
        angle = dot(v_obs_surf, vel_new) / (norm(v_obs_surf) * norm(vel_new))

        # project velocities along line of sight
        vs[i,j] = norm(vel_new) * angle

        # get surface element area
        dA = ρₛ^2.0 * sin(π/2.0 - ϕc[i]) * dϕ[i] * dθ[j]

        # get project surface element area
        dA_proj[i,j] = dA * abs(dot(v_obs_surf, xyz_new))
    end
end

# parse out coords
xyz_new = [M * i for i in get_xyz_for_surface.(ρₛ, ϕc, θc')]
xs_new = getindex.(xyz_new, 1)
ys_new = getindex.(xyz_new, 2)
zs_new = getindex.(xyz_new, 3)

# now plot the projection
idx = μs .< 0.0
μs_plot = copy(μs)
μs_plot[idx] .= NaN

vel_plot = copy(vs)
vel_plot[idx] .= NaN

# plt.pcolormesh(xs_new, zs_new, rad2deg.(acos.(μs_plot)))
# plt.pcolormesh(xs_new, zs_new, vs, cmap="seismic", vmin=-2000, vmax=2000, edgecolors="k")
plt.pcolormesh(xs_new, zs_new, dA_proj./sum(dA_proj))#, edgecolors="k")
plt.colorbar()
plt.show()

# cmap = mpl.cm.ScalarMappable(cmap="viridis")
# clrs = cmap.to_rgba(dA_proj)
# ax = plt.figure().add_subplot(projection="3d")
# ax.plot_surface(xs_new, ys_new, zs_new, facecolors=clrs)
# plt.xlabel("x")
# plt.ylabel("y")
# plt.zlabel("z")
# plt.show()
