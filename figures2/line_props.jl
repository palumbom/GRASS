# imports
using Pkg; Pkg.activate(".")
using CUDA
using GRASS
using Statistics

# plotting imports
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
using PyCall; animation = pyimport("matplotlib.animation");
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# get input data
lps = GRASS.LineProperties()
airwav = GRASS.get_rest_wavelength(lps)
llevel = GRASS.get_lower_level(lps)
ulevel = GRASS.get_upper_level(lps)
geff = GRASS.get_geff(lps)
height = GRASS.get_height(lps)
species = GRASS.get_species(lps)

cm = plt.cm.get_cmap("viridis")

fig = plt.figure()
ax1 = fig.add_subplot()
sc = ax1.scatter(airwav, ulevel, c=geff, cmap=cm)
for i in eachindex(airwav)
    if airwav[i] == 6302.4932
        continue
    end
    ax1.annotate(species[i], (airwav[i], ulevel[i] + 0.2), fontsize=10)
end
cb = fig.colorbar(sc)
cb.set_label(L"g_{\rm eff}")
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()
# ax1.set_xlim(xlim[1], xlim[2] + 50)
ax1.set_xlim(5200, 6400)
# ax1.set_ylim(ylim[1], ylim[2] + 0.3)
ax1.set_ylim(minimum(ulevel) - 0.5, maximum(ulevel) + 0.5)
ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
ax1.set_ylabel(L"{\rm Excitation\ Potential\ (eV)}")
fig.savefig(plotdir * "line_props.pdf")
plt.clf(); plt.close()
println(">>> Figure written to: " * plotdir * "line_props.pdf")