using GRASS
using SPICE
using CSV
using Statistics
using PyPlot
using DataFrames
using EchelleCCFs
using RvSpectMLBase
using EchelleInstruments

order_index = 61
full_pixels = 2048:7168
min_wav = 5435
max_wav = 5437

lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
airwav = lp.λrest[12]

path = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/neidL2_20231014T165429.fits"

function linear_interp(xs, ys; bc::T=NaN) where T<:Float64
    function f(x)
        if (((x < first(xs)) | (x > last(xs))) & !isnan(bc))
            return bc
        elseif x <= first(xs)
            return first(ys)
        elseif x >= last(xs)
            return last(ys)
        else
            i = searchsortedfirst(xs, x) - 1
            i0 = clamp(i, firstindex(ys), lastindex(ys))
            i1 = clamp(i+1, firstindex(ys), lastindex(ys))
            return (ys[i0] * (xs[i1] - x) + ys[i1] * (x - xs[i0])) / (xs[i1] - xs[i0])
        end
    end
    return f
end

function calc_line_bisector_at_abs_depth(λ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}; abs_depth::Real ) where { T1<:Real, T2<:Real }
    @assert length(λ) == length(flux)
    @assert 0.05 <= abs_depth <= 0.99
    idx_min_flux = argmin(flux)
    continuum = maximum(flux)
    if flux[idx_min_flux]/continuum > abs_depth  # Line isn't that deep!
        return nothing
    end
    target_flux = continuum*(1-abs_depth)
    idxhi = idx_min_flux-1+searchsortedfirst(view(flux,idx_min_flux:length(flux)), target_flux)
    idxlo = idxhi-1
    λ_hi = RvSpectMLBase.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
    idxlo = idx_min_flux+1-searchsortedfirst(view(flux,idx_min_flux:-1:1), target_flux )
    idxhi = idxlo+1
    λ_lo = RvSpectMLBase.interp_linear(x1=flux[idxlo],x2=flux[idxhi],y1=λ[idxlo],y2=λ[idxhi],xpred=target_flux)
    width = λ_hi-λ_lo
    return width
 end

#read in data
spectrum = NEID.read_data(path; normalization = :blaze)

#find where line is
chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav) 

chunk_flux = chunk.flux[pixels]

#interpolate line to be monotonic on each side
itp = linear_interp(chunk.λ[pixels], chunk_flux, bc=NaN)
flux_itp = itp.(LinRange(min_wav, max_wav, 100))

flux_itp .-= flux_itp[length(flux_itp)] 
flux_itp .+= 1 

#plot spectrum for specified order 
plt.figure()
plt.plot(spectrum.λ[full_pixels,order_index], spectrum.flux[full_pixels,order_index])
plt.plot(chunk.λ[pixels], chunk_flux)
plt.savefig("spectra.png")
plt.clf()

chunk_flux .-= chunk_flux[length(chunk_flux)] 
chunk_flux .+= 1 

#plot specific line of interest
plt.figure()
plt.plot(chunk.λ[pixels], chunk_flux)
plt.plot(LinRange(min_wav, max_wav, 100), flux_itp)
plt.axvline(x=airwav)
plt.savefig("spectra_zoom.png")

widths = calc_line_bisector_at_abs_depth(LinRange(min_wav, max_wav, 100), flux_itp, abs_depth = 0.3)

# --------------------------------------------------

# convolve IAG spectrum to LARS resolution
wavs_iag, flux_iag = GRASS.convolve_gauss(chunk.λ[pixels], chunk_flux, new_res=1e6, oversampling=4.0)
buff = 0.12575
idxl = findfirst(x -> x .>= airwav - buff, wavs_iag)
idxr = findfirst(x -> x .>= airwav + buff, wavs_iag)
iag_bot = minimum(view(flux_iag, idxl:idxr))
iag_depth = 1.0 - iag_bot

datadir = string(abspath("data"))
df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

# get the depth for the simulation
sim_depth = df[12, "optimized_depth"]

# simulate the spectrum
lines = [airwav]
depths = [sim_depth]
templates = ["FeI_5436.6"]
resolution = 7e5
spec = SpecParams(lines=lines, depths=depths, templates=templates,
                      resolution=resolution, buffer=1.5, oversampling=2.0)
N = 50

v_grid_cpu, ccf_cpu = GRASS.calc_ccf(LinRange(min_wav, max_wav, 100), flux_itp, spec)
rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)
println(rvs_cpu)
println(sigs_cpu)

#NEID location
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938 

neid_timestamps = ["2023-10-14T16:54:56.500000"]
#convert from utc to et as needed by SPICE
time_stamps = utc2et.(neid_timestamps)
Nt = length(time_stamps)
disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)

# simulate the spectrum
wavs_sim, flux_sim = GRASS.synthesize_spectra_eclipse(spec, disk, obs_long, obs_lat, alt, lines ./ 10.0, time_stamps, verbose=true, use_gpu=false)
flux_sim = dropdims(mean(flux_sim, dims=2), dims=2)

# interpolate iag on synth wavelength grid
itp = GRASS.linear_interp(wavs_iag, flux_iag)
flux_iag = itp.(wavs_sim)
wavs_iag = copy(wavs_sim)

    # get width in velocity for CCF
idxl_sim, idxr_sim = GRASS.find_wing_index(0.95, flux_sim)

    # get width in angstroms
width_ang = wavs_sim[idxr_sim] - wavs_sim[idxl_sim]

    # convert to velocity
width_vel = GRASS.c_ms * width_ang / wavs_sim[argmin(flux_sim)]
Δv_max = round((width_vel + 1e3)/100) * 100

v_grid_iag, ccf_iag = GRASS.calc_ccf(wavs_iag, flux_iag, lines, depths,
                                         7e5, Δv_step=100.0, Δv_max=Δv_max,
                                         mask_type=EchelleCCFs.GaussianCCFMask)

v_grid_sim, ccf_sim = GRASS.calc_ccf(wavs_sim, flux_sim, lines, depths,
                                         7e5, Δv_step=100.0, Δv_max=Δv_max,
                                         mask_type=EchelleCCFs.GaussianCCFMask)

# get bisectors
# vel_iag, int_iag = GRASS.calc_bisector(v_grid_iag, ccf_iag, nflux=50, top=top)
vel_sim, int_sim = GRASS.calc_bisector(v_grid_sim, ccf_sim, nflux=50, top=top)

# find mean velocities in order to align bisectors
N = 0.20
M = 0.70
idx1 = findfirst(x -> x .>= N * sim_depth + minimum(flux_sim), int_sim)
idx2 = findfirst(x -> x .>= M * sim_depth + minimum(flux_sim), int_sim)
if isnothing(idx2)
        idx2 = findfirst(x -> x .>= 0.9, int_sim)
end
rv_sim = mean(view(vel_sim, idx1:idx2))
println(rv_sim)

# idx1 = findfirst(x -> x .>= N * iag_depth + iag_bot, int_iag)
# idx2 = findfirst(x -> x .>= M * iag_depth + iag_bot, int_iag)
# if isnothing(idx2)
#         idx2 = findfirst(x -> x .>= 0.9, int_iag)
# end
# rv_iag = mean(view(vel_iag, idx1:idx2))
  
plt.figure()
plt.plot(wavs_sim, flux_sim, label = "sim")
plt.plot(wavs_iag, flux_iag)
plt.axvline(x=airwav)
plt.legend()
plt.savefig("test.png")                                        

v_grid_iag, ccf_iag = GRASS.calc_ccf(wavs_iag, flux_iag, spec) 
rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_iag, ccf_iag)

println(rvs_cpu)
println(sigs_cpu)