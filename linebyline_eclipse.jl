using GRASS
using SPICE
using CSV
using Statistics
using PyPlot
using DataFrames
using EchelleCCFs
using RvSpectMLBase
using EchelleInstruments

order_index = 75
full_pixels = 2048:7168
min_wav = 6169
max_wav = 6171

lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
airwav = lp.λrest[19]

path = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"

timestamps_full = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv", DataFrame)[!, "obsdate"]
timestamps = timestamps_full[16:length(timestamps_full)-150]

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

RV_list = Vector{Float64}(undef,length(timestamps)...)
for i in 1:length(timestamps)
    time = join(["neidL2_", timestamps[i][1], timestamps[i][2], timestamps[i][3], timestamps[i][4], timestamps[i][6], timestamps[i][7],
                timestamps[i][9], timestamps[i][10], "T", timestamps[i][12], timestamps[i][13], timestamps[i][15], timestamps[i][16],
                timestamps[i][18], timestamps[i][19], ".fits"])

    spectrum = NEID.read_data(joinpath(path, time); normalization = :blaze)

    #find where line is
    chunk = NEID.ChunkOfSpectrum(spectrum, order_index, full_pixels) 
    pixels = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav) 

    chunk_flux = chunk.flux[pixels]

    #interpolate line to be monotonic on each side
    itp = linear_interp(chunk.λ[pixels], chunk_flux, bc=NaN)
    flux_itp = itp.(LinRange(min_wav, max_wav, 100))

    #plot spectrum for specified order 
    plt.figure()
    plt.plot(spectrum.λ[full_pixels,order_index], spectrum.flux[full_pixels,order_index])
    plt.plot(chunk.λ[pixels], chunk_flux)
    plt.savefig("spectra.png")
    plt.clf()

    #plot specific line of interest
    plt.figure()
    plt.plot(chunk.λ[pixels], chunk_flux)
    plt.plot(LinRange(min_wav, max_wav, 100), flux_itp)
    plt.axvline(x=airwav)
    plt.savefig("spectra_zoom.png")

    datadir = string(abspath("data"))
    df = CSV.read(joinpath(datadir, "optimized_depth.csv"), DataFrame)

    # get the depth for the simulation
    sim_depth = df[19, "optimized_depth"]

    # simulate the spectrum
    lines = [airwav]
    depths = [sim_depth]
    templates = ["FeI_5436.6"]
    resolution = 7e5
    spec = SpecParams(lines=lines, depths=depths, templates=templates,
                        resolution=resolution, buffer=1.5, oversampling=2.0)

    v_grid_cpu, ccf_cpu = GRASS.calc_ccf(LinRange(min_wav, max_wav, 100), flux_itp, spec)
    rvs_cpu, sigs_cpu = GRASS.calc_rvs_from_ccf(v_grid_cpu, ccf_cpu)

    RV_list[i] = rvs_cpu
end

# RV_list .-= RV_list[length(RV_list)] 
# RV_list .+= 1 

println(RV_list)

plt.figure()
plt.plot(1:length(timestamps), RV_list)
plt.savefig("rm.png")