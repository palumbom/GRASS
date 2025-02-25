using Revise
using CSV
using CUDA
using Revise, GRASS
using SPICE
using Dates
using DataFrames

# plotting impors
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
using LaTeXStrings

# get python imports 
wcs = pyimport("astropy.wcs")
SkyCoord = pyimport("astropy.coordinates").SkyCoord
RegionGeom = pyimport("gammapy.maps").RegionGeom
CircleSkyRegion = pyimport("regions").CircleSkyRegion
WCSAxes = pyimport("astropy.visualization.wcsaxes").WCSAxes
coord = pyimport("astropy.coordinates")
u = pyimport("astropy.units")
AltAz = pyimport("astropy.coordinates").AltAz
EarthLocation = pyimport("astropy.coordinates").EarthLocation
apTime = pyimport("astropy.time").Time

# get kernels
GRASS.get_kernels()

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

# get time_stamps CHANGE
time_stamps = range(utc2et("2026-01-10T07:00:00.0"), utc2et.("2026-01-10T12:00:00.0"), step = 105.0)
# time_stamps = range(utc2et("2023-10-14T14:00:00.0"), utc2et.("2023-10-14T19:00:00.0"), step = 100.0)
# time_stamps_utc = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
# time_stamps = utc2et.(time_stamps_utc)

# set up disk params
N = 300
Nt = length(time_stamps)
Nsubgrid = 10
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
LD_type = "SSD"
extinction_coeff = DataFrame(CSV.File(joinpath(GRASS.datdir, "NEID_three_ext_coeff_KSSD.csv")))
neid_ext_coeff = extinction_coeff[extinction_coeff[!, "Wavelength"] .== lines[1], "Ext1"]

# allocate memory for gpu
gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, length(spec.lines))

# create a subdir to write to 
outdir = joinpath(pwd(), "icymoon_frames") #CHANGE
if !isdir(outdir)
    mkdir(outdir)
end

# loop over timestamps 
for t in eachindex(time_stamps)
    # get the geometry data
    GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, obs_alt, lines, 
                                       LD_type, 0.0, neid_ext_coeff[1], disk, gpu_allocs) #CHANGE

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
    dat[isnan.(μs)] .= NaN

    # create a WCS for the projection 
    # w = wcs.WCS(naxis=2)
    # w.wcs.ctype = ["RA---ZEA", "DEC--ZEA"]

    # get the earth location 
    # earth_loc = EarthLocation(lat=obs_lat * u.deg, lon=obs_long * u.deg, height = obs_alt * u.km)

    the_coord = SkyCoord(az[1], alt[1], unit="deg", frame="icrs")
    # the_frame = AltAz(obstime=apTime(time_stamps_utc[1]), location=earth_loc)
    the_wcs = RegionGeom(CircleSkyRegion(center=the_coord, radius=0.5 * u.deg)).wcs

    # the_wcs = wcs.WCS(naxis=2)
    # the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]  # Azimuth and Altitude
    # the_wcs.wcs.crval = [180, 90]  # Reference point
    # the_wcs.wcs.cdelt = [-1, 1]  # Scaling (degrees/pixel)
    # the_wcs.wcs.crpix = [180, 45]  # Reference pixel

    # set up figure
    fig = plt.figure()
    ax1 = plt.subplot(projection=the_wcs)

    # set unit
    # ax1.coords[1].set_format_unit("deg", show_decimal_unit=true)
    # ax1.coords[2].set_format_unit("deg", show_decimal_unit=true)
    ax1.coords[1].set_major_formatter("d.dd")
    ax1.coords[2].set_major_formatter("d.dd")    
    # ax1.coords[1].set_axislabel("Azimuth (deg)")
    # ax1.coords[2].set_axislabel("Altitude (deg)")
    ax1.set_xlabel("Azimuth (deg)")
    ax1.set_ylabel("Altitude (deg)")

    # make colormap for data
    cmap = plt.cm.seismic
    # norm = mpl.colors.Normalize(vmin=minimum(dat), vmax=maximum(dat))
    norm = mpl.colors.Normalize(vmin=-2400, vmax=2400)
    smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    # add a colorbar 
    fig.colorbar(smap, ax=ax1, orientation="vertical", label="Velocity (m/s)")

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
    fig.savefig(joinpath(outdir, "eclipse_" * string(t) * ".png"))
    plt.close()
end

### then use gifski on command line to stitch together pngs
### download the cli exe here: https://gif.ski
# gifski -o output.gif eclipse_*.png --fps=8