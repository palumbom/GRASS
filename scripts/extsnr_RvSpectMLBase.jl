using CSV
using JLD2
using GRASS
using SPICE
using FITSIO
using FileIO
using Statistics
using DataFrames
using EchelleCCFs
using RvSpectMLBase
using EchelleInstruments
using Interpolations
using EchelleCCFs: λ_air_to_vac, λ_vac_to_air
using Dierckx
using LsqFit
# plotting
using LaTeXStrings
import PyPlot
plt = PyPlot
mpl = plt.matplotlib
using Colors, ColorSchemes
colors = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]

using PyCall
py"""
import sys
sys.path.append('.')
"""
np = pyimport("numpy")
pd = pyimport("pandas")
scipy_interp = pyimport("scipy.interpolate")
CubicSpline = scipy_interp.CubicSpline
mdates = pyimport("matplotlib.dates")

# determine lines to be considered
lp = GRASS.LineProperties(exclude=["CI_5380", "NaI_5896"])
line_names = GRASS.get_name(lp)
airwav = lp.λrest
vacwav = λ_air_to_vac.(airwav)

cmap = ColorSchemes.coolwarm   
normalize01(x) = (maximum(vacwav) == minimum(vacwav) ? 0.5 : clamp((x - minimum(vacwav))/(maximum(vacwav) - minimum(vacwav)), 0.0, 1.0))

# line information for identification
orders = [57, 57, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 64, 64, 74, 74, 75, 75, 75, 75, 77, 77]
airmass = [2.5632073517047638, 2.5361097748885557, 2.5093661650400723, 2.4832872029346307, 2.458152719768879, 2.4333280353218156, 2.4091026139254765, 2.3857380385003344, 2.362645575113234, 2.3400953849743926, 2.3183323449274487, 2.2968091303568645, 2.27577829301206, 2.2554693819371017, 2.2353725338593415, 2.2157241944955177, 2.1967398191699803, 2.1779435542209677, 2.159557077563393, 2.1415684495759097, 2.1241760200750517, 2.106944668512874, 2.0900782010965315, 2.0737632866181173, 2.05759228654005, 2.0417568347767414, 2.0264327849819934, 2.011237707679954, 1.9963519894425743, 1.9819414507550772, 1.9676468497988242, 1.9536381582317257, 1.9400718707467606, 1.926610129678612, 1.9134132425208898, 1.9006290194999498, 1.8879393729586693, 1.8754956248574794, 1.8634374833021545, 1.8514651818936574, 1.8397216717571738, 1.8282016018286706, 1.817034683533013, 1.8059435946574782, 1.7950609742097614, 1.7845095958663924, 1.7740276105655446, 1.763740511098882, 1.753764560828184, 1.743852354616247, 1.7341226871090905, 1.7246857022512947, 1.715307505275213, 1.7061006017579452, 1.697169374222985, 1.688292589769291, 1.6795768389575176, 1.6711210712224671, 1.6627159428484009, 1.6544624719613135, 1.646454476321036, 1.638493796560436, 1.630676192498115, 1.6230906273882857, 1.6155494821322123, 1.6081435466952902, 1.6008703361742513, 1.5938127284826231, 1.5867962627606653, 1.5799055052643187, 1.573219043936084, 1.5665716579824998, 1.5600435303724371, 1.5537091700100665, 1.5474121144362678, 1.5412283821572987, 1.5352285977595637, 1.5292646143956121, 1.523408487955313, 1.5177271385876636, 1.5120803277757258, 1.5065363351237604, 1.5011585426132867, 1.4958142450515572, 1.4905681178810164, 1.4854801586717268, 1.4804248500221022, 1.4754634222942877, 1.4705945147488082, 1.4658738276317298, 1.4611849448608079, 1.456584725994602, 1.4521258025734258, 1.447698240609205, 1.4433557809843884, 1.43914810522267, 1.43497149743696, 1.4308767024770255, 1.4269105616726943, 1.4229753336016377, 1.419118880539379, 1.4153853052103214, 1.4116826158282119, 1.408055897033766, 1.4045466068708823, 1.4010682957774412, 1.3976633677805828, 1.3943707230207358, 1.3911092629152901, 1.387918800523127, 1.3847985802902372, 1.381784212071226, 1.3788014670576347, 1.3758868280098038, 1.3730735241834369, 1.3702922880917996, 1.3675771947689528, 1.364959158943391, 1.3623737296001717, 1.3598526423845454, 1.3574245581270798, 1.3550297095145418, 1.3526975555563139, 1.350454558686985, 1.3482455138644898, 1.346097660961062, 1.3440353146548754, 1.34200772145745, 1.3400399546266697, 1.3381542271029911, 1.3363041363785326, 1.334512636897446, 1.3327998814765274, 1.331123727652466, 1.329505054361501, 1.3279435167148383, 1.3264565754953044, 1.3250076481928261, 1.323614903675427, 1.3222938236556199, 1.3210119084301528, 1.319785336067303, 1.3186276374092976, 1.3175103337752088, 1.3164476420975404, 1.3154511668285525, 1.3144963965600829, 1.3135956129740973, 1.3127585150081922, 1.3119645125498263, 1.3112239742717597, 1.310544710948751, 1.3099100153508805, 1.3093283618247469, 1.3088056868926414, 1.3083291351789732, 1.307905301914301, 1.3075341005197691]
airmass = airmass[1:130]
full_pixels = 2048:7168 

grass_intensity_array = nothing
jldopen("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/SSD/data/projected_SSD_4parameter_gpu.jld2", "r") do file
    global grass_intensity_array
    grass_intensity_array = (file["intensity"])
end

path = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
timestamps_full_october = CSV.read("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv", DataFrame)[!, "filename"]
timestamps = timestamps_full_october[16:length(timestamps_full_october)-178]

neid_intensity_array = Vector{Vector{Float64}}(undef, length(line_names))
for i in 1:length(line_names)
    time_array = Vector{Float64}(undef, length(timestamps))
    for j in 1:length(timestamps)
        order = orders[i]
        min_wav = vacwav[i] - 5.0
        max_wav = vacwav[i] + 5.0
        buffer = 0.2

        spectrum = NEID.read_data(joinpath(path, timestamps[j]))
        chunk = NEID.ChunkOfSpectrum(spectrum, order, full_pixels) 
        pixel_normal = NEID.find_pixels_for_line_in_chunk(chunk, min_wav, max_wav)
        pixels_inner = NEID.find_pixels_for_line_in_chunk(chunk, vacwav[i] - buffer, vacwav[i] + buffer)
        chunk_flux = chunk.flux[pixels_inner]
        time_array[j] = np.nanmedian(chunk_flux)
    end
    neid_intensity_array[i] = time_array
end

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for (i, value) in enumerate(vacwav)
#     color = get(ColorSchemes.coolwarm, normalize01(value)) 
#     color = (red(color), green(color), blue(color), 1.0)
#     ax1.scatter(airmass, grass_intensity_array[i] ./ maximum(grass_intensity_array[i]), color=color, s = 1)
#     ax1.plot(airmass, neid_intensity_array[i] ./ maximum(neid_intensity_array[i]), color=color)
# end
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux/intensity") 
# ax1.axvline(x = airmass[49])
# plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/investigations/extinction/no_ext/flux_comp_SSD_new.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for (i, value) in enumerate(vacwav)
#     color = get(ColorSchemes.coolwarm, normalize01(value)) 
#     color = (red(color), green(color), blue(color), 1.0)
#     grass_normalized = grass_intensity_array[i] ./ maximum(grass_intensity_array[i])
#     neid_normalized = neid_intensity_array[i] ./ maximum(neid_intensity_array[i])
#     comp_arr = neid_normalized ./ grass_normalized
#     ax1.plot(airmass, comp_arr, color=color)
# end
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux/intensity") 
# ax1.axvline(x = airmass[49])
# plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/investigations/extinction/no_ext//coeff_eclipse_SSD_new.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()  
# for (i, value) in enumerate(vacwav)
#     color = get(ColorSchemes.coolwarm, normalize01(value)) 
#     color = (red(color), green(color), blue(color), 1.0)
#     grass_normalized = grass_intensity_array[i] ./ maximum(grass_intensity_array[i])
#     neid_normalized = neid_intensity_array[i] ./ maximum(neid_intensity_array[i])
#     comp_arr = neid_normalized - grass_normalized
#     ax1.scatter(airmass, comp_arr, color=color, s = 1)  
# end
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux") 
# ax1.axvline(x = airmass[49])
# plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/investigations/extinction/no_ext//flux_comp_diff_SSD_new.png")
# plt.clf()

fig = plt.figure()
ax1 = fig.add_subplot()
p0 = [0.1]
ext_coeff = Vector{Float64}(undef, length(line_names))
for (i, value) in enumerate(vacwav)
    model(x, p) = grass_intensity_array[i] .* np.exp(-x .* p[1]) ./ maximum(grass_intensity_array[i] .* np.exp(-x .* p[1]))
    fit_model = curve_fit(model, airmass, neid_intensity_array[i] ./ maximum(neid_intensity_array[i]), p0)

    color = get(ColorSchemes.coolwarm, normalize01(value)) 
    color = (red(color), green(color), blue(color), 1.0)
    new_grass_model = model(airmass, fit_model.param[1])
    ax1.scatter(airmass, new_grass_model, color=color, s = 1)
    ax1.plot(airmass, neid_intensity_array[i] ./ maximum(neid_intensity_array[i]), color=color)

    ext_coeff[i] = fit_model.param[1]
end
println(ext_coeff)
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux/intensity") 
# ax1.axvline(x = airmass[49])
# plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/investigations/extinction/ext/flux_comp_SSD_new.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for (i, value) in enumerate(vacwav)
#     color = get(ColorSchemes.coolwarm, normalize01(value)) 
#     color = (red(color), green(color), blue(color), 1.0)
#     model(x, p) = grass_intensity_array[i] .* np.exp(-x .* p[1]) ./ maximum(grass_intensity_array[i] .* np.exp(-x .* p[1]))
#     grass_normalized = model(airmass, ext_coeff[1])
#     neid_normalized = neid_intensity_array[i] ./ maximum(neid_intensity_array[i])
#     comp_arr = neid_normalized ./ grass_normalized
#     ax1.plot(airmass, comp_arr, color=color)
# end
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux/intensity") 
# ax1.axvline(x = airmass[49])
# plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/investigations/extinction/ext/coeff_eclipse_SSD_new.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()  
# for (i, value) in enumerate(vacwav)
#     color = get(ColorSchemes.coolwarm, normalize01(value)) 
#     color = (red(color), green(color), blue(color), 1.0)
#     model(x, p) = grass_intensity_array[i] .* np.exp(-x .* p[1]) ./ maximum(grass_intensity_array[i] .* np.exp(-x .* p[1]))
#     grass_normalized = model(airmass, ext_coeff[1])
#     neid_normalized = neid_intensity_array[i] ./ maximum(neid_intensity_array[i])
#     comp_arr = neid_normalized - grass_normalized
#     ax1.scatter(airmass, comp_arr, color=color, s = 1)  
# end
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux") 
# ax1.axvline(x = airmass[49])
# plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/investigations/extinction/ext/flux_comp_diff_SSD_new.png")
# plt.clf()