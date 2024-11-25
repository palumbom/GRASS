using Revise
using CSV
using CUDA
using GRASS
using SPICE
using Dates
using DataFrames

# plotting impors
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
using LaTeXStrings

# get python imports 
# wcs = pyimport("astropy.wcs")
# SkyCoord = pyimport("astropy.coordinates").SkyCoord
# AltAz = pyimport("astropy.coordinates").AltAz
# EarthLocation = pyimport("astropy.coordinates").EarthLocation

function get_gmst(jd_ut1)
    t = (jd_ut1 - 2451545.0) / 36525.0
    era = earth_rotation_angle(jd_ut1)

    gmst = (era + (0.014506 + 4612.15739966*t + 1.39667721*t^2.0 - 0.00009344*t^3.0 + 0.00001882*t^4.0) / 60.0 / 60.0 * π / 180.0) % (2π)
    if gmst < 0.0
        gmst += 2π
    end

    return gmst
end

function earth_rotation_angle(jd_ut1)
    t = jd_ut1 - 2451545.0;
    frac = jd_ut1 % 1.0;
    
    era = (2π * (0.7790572732640 + 0.00273781191135448 * t + frac)) % (2π);
    if era < 0.0
        era += 2π
    end
    
    return era
end

function get_alt_az(et, ra, dec, obs_lat, obs_long)
    # get greenwhich mean sidereal time 
    jd = tryparse(Float64, et2utc(et, "J", 10)[4:end])
    gmst = rad2deg(get_gmst(jd))

    # get local sidereal time 
    lst = mod(gmst + obs_long, 360.0) 

    # get hour angle
    ha = mod(lst - ra, 360.0)

    # quantities in rad 
    ha_rad = deg2rad(ha)
    dec_rad = deg2rad(dec)
    lat_rad = deg2rad(obs_lat)

    # alt and az
    alt = asin(sin(lat_rad) * sin(dec_rad) + cos(lat_rad) * cos(dec_rad) * cos(ha_rad))

    num = -sin(ha_rad) 
    den = tan(dec_rad) * cos(lat_rad) - sin(lat_rad) * cos(ha_rad)
    azimuth = atan(num, den)
    if azimuth < 0.0
        azimuth += 2π
    end

    # convert to deg and return
    return rad2deg(alt), rad2deg(azimuth)
end

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
obs_alt = 2.097938

# get time_stamps
time_stamps = range(utc2et("2023-10-14T14:00:00.0"), utc2et.("2023-10-14T19:00:00.0"), step = 100.0)

# set up disk params
N = 197
Nt = length(time_stamps)
Nsubgrid = 40
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=Nsubgrid)

# stuff for line
lines = [5434.5232] # array of line centers 
depths = [0.6]   # array of line depths
templates = ["FeI_5434"] # template data to use
variability = trues(length(lines))  # whether or not the bisectors should "dance"
blueshifts = zeros(length(lines))   # set convective blueshift value
resolution = 7e5                    # spectral resolution

# make the disk and spec composite type instances
spec = GRASS.SpecParams(lines=lines, depths=depths, variability=variability,
                        blueshifts=blueshifts, templates=templates, 
                        resolution=resolution) 

# limb darkening and extinction
LD_type = "KSSD"
extinction_coeff = DataFrame(CSV.File(joinpath(GRASS.datdir, "NEID_three_ext_coeff_KSSD.csv")))
neid_ext_coeff = extinction_coeff[extinction_coeff[!, "Wavelength"] .== lines[1], "Ext1"]

# allocate memory for gpu
gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, length(spec.lines))

# create a WCS for the projection 
# w = wcs.WCS(naxis=2)
# w.wcs.ctype = ["RA---ZEA", "DEC--ZEA"]

# initialize figure objects outside time loop
fig = plt.figure()
ax1 = plt.subplot()#projection=w)

# make colormap for data
cmap = plt.cm.seismic
# norm = mpl.colors.Normalize(vmin=minimum(dat), vmax=maximum(dat))
norm = mpl.colors.Normalize(vmin=-2400, vmax=2400)
smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

# add a colorbar 
fig.colorbar(smap, ax=ax1, orientation="vertical", label="Velocity (m/s)")

# set labels
ax1.set_xlabel(L"{\rm Azimuth\ (deg)}", labelpad=10)
ax1.set_ylabel(L"{\rm Altitide\ (deg)}", labelpad=10)

# set limits
# ax1.set_xlim(100.0, 210.0)
# ax1.set_ylim(0.0, 50.0)
# ax1.set_aspect("equal")

# ax1.set_xlim(reverse(ax1.get_xlim())...)
# ax1.set_xlabel(L"{\rm Right\ Ascension}")
# ax1.set_ylabel(L"{\rm Declination}")

# function to loop over time
function sun_at_time(t)
# t = length(time_stamps) - 1
    # python is zero-indexed
    t += 1

    # get the geometry data
    GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, obs_alt, lines, 
                                       LD_type, 1.0, neid_ext_coeff[1], disk, gpu_allocs)

    # parse out data
    μs = Array(gpu_allocs.μs)
    dA = Array(gpu_allocs.dA)
    ld = Array(gpu_allocs.ld)[:,:,1]
    vrot = Array(gpu_allocs.z_rot) .* GRASS.c_ms
    ext = Array(gpu_allocs.ext)[:,:,1]

    # get grid edges for plot
    xs, ys, zs = GRASS.calc_grid_edge_xyz(time_stamps[t], obs_long, obs_lat, obs_alt, disk, gpu_allocs)

    # get ra and dec
    out = SPICE.recrad.(zip(xs, ys, zs))
    rr = getindex.(out, 1)
    ra = rad2deg.(getindex.(out, 2))
    dec = rad2deg.(getindex.(out, 3))

    # get alt and az 
    out = get_alt_az.(time_stamps[t], ra, dec, obs_lat, obs_long)
    alt = getindex.(out, 1)
    az = getindex.(out, 2)

    # mask coords 
    idx = rr .<= 0.0
    ra[idx] .= NaN
    dec[idx] .= NaN
    alt[idx] .= NaN
    az[idx] .= NaN

    # select data
    dat = deepcopy(vrot)

    # mask data
    dat[μs .<= 0.0] .= NaN

    # clear collections if they exist
    if t > 1 
        for artist in ax1.patches
            artist.remove()
        end
    end

    aa = zeros(4)
    bb = zeros(4)

    cell_order = [1,2,4,3]

    # loop over and plot grid positions
    fills = []
    for i in 1:size(dat, 1)
        for j in 1:size(dat, 2)
            # get m and n 
            m = 2 * (i - 1) + 1
            n = j

            if isnan(dat[i,j])
                continue
            end

            # check that we haven't overflowed
            if any(iszero.(ra[m:m+1,n:n+1]))
                continue
            end

            aa .= vec(az[m:m+1,n:n+1])
            bb .= vec(alt[m:m+1,n:n+1])

            aa .= view(aa, cell_order)
            bb .= view(bb, cell_order)

            this_fill = ax1.fill(aa, bb, c=smap.to_rgba(dat[i,j]))
            push!(fills, this_fill)
        end
    end

    # get limits
    buff = 0.05
    ax1.set_xlim(minimum(az[.!idx]) - buff, maximum(az[.!idx]) + buff)
    ax1.set_ylim(minimum(alt[.!idx]) - buff, maximum(alt[.!idx]) + buff)
    ax1.set_aspect("equal")#, "datalim")
    return nothing
end

ani = animation.FuncAnimation(fig, sun_at_time, frames=length(time_stamps), interval=1)
ani.save(filename="eclipse2.gif", writer="pillow", fps=8)
ani.event_source.stop()
plt.close("all")

# plt.show()
