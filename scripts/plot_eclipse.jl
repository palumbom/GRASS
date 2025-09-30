using Revise
using CSV
using CUDA
using Revise, GRASS
using SPICE
using Dates
using DataFrames
using StaticArrays
using JLD2

# plotting impors
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
using LaTeXStrings

# Load sunspots from detection
sunspot_coords = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/investigations/sunspots/sunspots.csv", DataFrame)

const Rsun_km = 695700.0    # Solar radius (km). Replace with SPICE.bodvrd if you want kernel-derived value.

# Convert Carrington lon/lat (deg) -> RA/Dec (deg) as seen from Earth center at ET `et`.
function carrington_to_radec_as_seen_from_earth(lon_deg, lat_deg, et)
    # 1) body-fixed vector in IAU_SUN frame (km)
    lon = deg2rad(lon_deg)
    lat = deg2rad(lat_deg)
    r_body = SVector{3}(Rsun_km * cos(lat) * cos(lon),
                       Rsun_km * cos(lat) * sin(lon),
                       Rsun_km * sin(lat))

    # 2) rotate to J2000 at epoch 'et'
    # pxform returns rotation matrix R such that: r_j2000 = R * r_body
    R = SPICE.pxform("IAU_SUN", "J2000", et)
    r_j2000 = R * r_body

    # 3) get Earth's position w.r.t. Sun (vector from Sun -> Earth in J2000)
    # spkpos returns (pos, lt) in many SPICE bindings; grab the first element (position vector)
    spk = SPICE.spkpos("EARTH", et, "J2000", "NONE", "SUN")
    r_earth = getindex(spk, 1)   # should be a length-3 array
    r_earth = SVector{3}(r_earth...)   # make immutable vector for arithmetic

    # 4) vector from Earth -> spot = r_spot_j2000 - r_earth_j2000
    r_e2s = r_j2000 .- r_earth

    # 5) spherical coords (radius, RA, Dec) of that vector
    out = SPICE.recrad(r_e2s)   # returns (r, lon, lat) in radians where lon ~ RA, lat ~ Dec
    ra_deg  = rad2deg(getindex(out, 2))
    dec_deg = rad2deg(getindex(out, 3))

    return ra_deg, dec_deg, getindex(out,1)  # return radius (for visibility checks) too
end

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
time_stamps_utc = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
time_stamps = utc2et.(time_stamps_utc[1:length(time_stamps_utc)-28])

solar_horizon_angle = [13.69214256, 13.83337338, 13.97779268, 14.12369855, 14.26931778, 14.41821379,
 14.56863172, 14.71874374, 14.87222317, 15.02726097, 15.18197332, 15.34014609,
 15.49991472, 15.65933803, 15.82231716, 15.98693064, 16.15117858, 16.31908009,
 16.48865543, 16.65991802, 16.83078732, 17.00544436, 17.18182936, 17.35779963,
 17.53766089, 17.71929168, 17.90048562, 18.0856761,  18.27267846, 18.45922109,
 18.64986805, 18.84236995, 19.03438833, 19.23062102, 19.42875231, 19.62637529,
 19.82832467, 20.03221677, 20.23557463, 20.44337295, 20.65315849, 20.8649462,
 21.07616302, 21.29197539, 21.50983468, 21.72709386, 21.94906676, 22.17313129,
 22.39656481, 22.62483167, 22.85523465, 23.08497384, 23.31966724, 23.5565408,
 23.79271569, 24.03396663, 24.27744111, 24.5201797,  24.7681169,  25.01832006,
 25.26774752, 25.52249654, 25.77955268, 26.0357905,  26.29747282, 26.56150185,
 26.82789041, 27.09339885, 27.36451523, 27.63802805, 27.91061039, 28.18892203,
 28.46966454, 28.74942238, 29.03502939, 29.32309881, 29.61012523, 29.90311863,
 30.19860264, 30.4929808,  30.79344113, 31.09641642, 31.39821816, 31.70621399,
 32.01674479, 32.32602911, 32.64161548, 32.95975195, 33.28044239, 33.59978027,
 33.92555687, 34.25389518, 34.58079442, 34.91422823, 35.25022529, 35.58469034,
 35.92577932, 36.26942629, 36.61144158, 36.96016284, 37.31142946, 37.66095766,
 38.0172657,  38.37609853, 38.73307886, 39.09690389, 39.46322465, 39.82757125,
 40.19881757, 40.57252151, 40.94866823, 41.32266629, 41.70362072, 42.08696704,
 42.46802471, 42.85606762, 43.24644145, 43.6343787,  44.02931761, 44.42651603,
 44.82112206, 45.22273306, 45.62652147, 46.02755394, 46.43558098, 46.8456925,
 47.25287707, 47.66703162, 48.0831668,  48.49619692]

# set up disk params
N = 40
Nt = length(time_stamps)
Nsubgrid = 2
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
LD_type = "SSD_4parameter"
extinction_coeff = DataFrame(CSV.File(joinpath(GRASS.datdir, "NEID_two_ext_SSD_4parameter.csv")))

# allocate memory for gpu
gpu_allocs = GRASS.GPUAllocsEclipse(spec, disk, length(spec.lines))

# create a subdir to write to 
outdir = joinpath(pwd(), "eclipse_plots")
if !isdir(outdir)
    mkdir(outdir)
end

# loop over timestamps 
spots_ra_time_arr = Vector{Vector{Float64}}(undef,size(time_stamps)...)
spots_dec_time_arr = Vector{Vector{Float64}}(undef,size(time_stamps)...)
for t in eachindex(time_stamps)
    if t < 48
        neid_ext_coeff = extinction_coeff[extinction_coeff[!, "Wavelength"] .== lines[1], "ext1"][1]
    elseif t >= 48 
        neid_ext_coeff = extinction_coeff[extinction_coeff[!, "Wavelength"] .== lines[1], "ext2"][1]
    end
    # get the geometry data
    GRASS.calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, obs_alt, lines, 
                                       LD_type, 1.0, neid_ext_coeff, disk, gpu_allocs)

    # parse out data
    μs = Array(gpu_allocs.μs)
    dA = Array(gpu_allocs.dA)
    ld = Array(gpu_allocs.ld)[:,:,1]
    vrot = Array(gpu_allocs.z_rot) .* GRASS.c_ms
    ext = Array(gpu_allocs.ext)[:,:,1]

    # vrot = Array(gpu_allocs.earth_v)
    # println(minimum(vrot))
    # println(maximum(filter(x -> x < 0, vrot)))
    # println(abs(minimum(vrot))-abs(maximum(filter(x -> x < 0, vrot))))

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

    the_coord = SkyCoord(az[1], alt[1], unit="deg", frame="icrs")
    the_wcs = RegionGeom(CircleSkyRegion(center=the_coord, radius=0.5 * u.deg)).wcs

    # set up figure
    fig = plt.figure(figsize=(8, 8))
    ax1 = plt.subplot(projection=the_wcs)

    ax1.coords[1].set_major_formatter("d.dd")
    ax1.coords[2].set_major_formatter("d.dd")   
    ax1.tick_params(axis="both", labelsize=14) 
    ax1.set_xlabel("Azimuth", fontsize = 14)
    ax1.set_ylabel("Altitude", fontsize = 14)

    # make colormap for data
    cmap = plt.cm.seismic
    norm = mpl.colors.Normalize(vmin=-2400, vmax=2400)
    smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    # add a colorbar 
    cbar = fig.colorbar(smap, ax=ax1, orientation="vertical")
    cbar.set_label(label="Velocity (m/s)", size = 14)
    cbar.ax.tick_params(labelsize=14)

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
    split_time = split(split(time_stamps_utc[t], "T")[2], ":")
    ax1.text(x=0.02, y=0.98, s=split_time[1] * ":" * split_time[2], transform=ax1.transAxes,
        horizontalalignment="left", verticalalignment="top", fontsize=14)
    ax1.text(x=0.98, y=0.98, s="$(round(solar_horizon_angle[t], digits = 2))°", transform=ax1.transAxes,
        horizontalalignment="right", verticalalignment="top", fontsize=14)
    # plax1.set_aspect('equal')

    spot_radec = Vector{Tuple{Float64,Float64,Float64}}()
    spots_ra_arr = Vector{Float64}(undef,size(sunspot_coords,1)...)
    spots_dec_arr = Vector{Float64}(undef,size(sunspot_coords,1)...)
    for (i, row) in enumerate(eachrow(sunspot_coords))
        lon = row[:lon]; lat = row[:lat]
        ra, dec, rr = carrington_to_radec_as_seen_from_earth(lon, lat, time_stamps[t])
        spots_ra_arr[i] = ra
        spots_dec_arr[i] = dec
        push!(spot_radec, (ra, dec, rr))
    end
    spots_ra_time_arr[t] = spots_ra_arr
    spots_dec_time_arr[t] = spots_dec_arr

    # (2) convert RA/Dec -> Alt/Az and plot markers (note: your plotting coord system uses az,alt)
    for (ra, dec, rr) in spot_radec
        # optional: skip spots on far side / behind Earth by checking rr or geometry
        # if rr <= 0.0 continue end   # usually rr>0 always; better visibility check follows

        alt, az = get_alt_az(time_stamps[t], ra, dec, obs_lat, obs_long)
        # optional: skip NaNs or obviously below horizon:
        if !isnan(alt) && !isnan(az)
            # you use az as X, alt as Y when drawing grid cells, so do the same:
            ax1.plot(az, alt, "x", markersize=8, markeredgewidth=2, color = "k")
        end
    end

    fig.savefig(joinpath(outdir, "eclipse_" * string(t) * ".png"), dpi = 100, bbox_inches="tight")
    plt.close()
end

@save "data/sunspots_ra_dec.jld2"
jldopen("data/sunspots_ra_dec.jld2", "a+") do file
    file["spots_ra_time_arr"] = deepcopy(spots_ra_time_arr) 
    file["spots_dec_time_arr"] = deepcopy(spots_dec_time_arr) 
end