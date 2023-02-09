using Pkg; Pkg.activate(".")
using CSV
using CUDA
using GRASS
using Printf
using LsqFit
using DataFrames
using Statistics
using EchelleCCFs
using Distributions

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "figures1/fig.mplstyle")

# get command line args and output directories
run, plot = parse_args(ARGS)
grassdir, plotdir, datadir = check_plot_dirs()

# determine whether to use GPU
use_gpu = CUDA.functional()

function blueshift_vs_depth(depths::AbstractArray{Float64,1};
                            use_gpu::Bool=use_gpu,
                            blueshifts=[])
    # set up stuff for lines
    lines = [5434.5]
    templates = ["FeI_5434"]
    resolution = 7e5
    disk = DiskParams(N=132, Nt=50)

    # allocate memory and loop over depths
    rvs_avg = zeros(length(depths))
    std_avg = zeros(length(depths))
    for (idx, depth) in enumerate(depths)
        @printf(">>> Doing depth %.2f \r", depth)
        if isempty(blueshifts)
            spec = SpecParams(lines=lines, depths=[depth], templates=templates, resolution=resolution)
        else
            spec = SpecParams(lines=lines, depths=[depth], blueshifts=[blueshifts[idx]], templates=templates, resolution=resolution)
        end

        # synthesize spectra
        lambdas1, outspec1 = synthesize_spectra(spec, disk, use_gpu=use_gpu, verbose=false)

        # calculate the RVs
        v_grid, ccf1 = calc_ccf(lambdas1, outspec1, spec)
        rvs1, sigs1 = calc_rvs_from_ccf(v_grid, ccf1)
        rvs_avg[idx] = mean(rvs1)
        std_avg[idx] = std(rvs1)
    end
    println()
    return rvs_avg, std_avg
end

# read in the data table from IAG paper
df_iag = GRASS.read_iag_blueshifts()
waves = df_iag.wavelength
depths = df_iag.depth
blueshift = df_iag.blueshift

# get indices for bins
bin_edges = range(0.025, 1.0-0.025, step=0.05)
bin_centers = bin_edges[1:end-1] .+ 0.5 * step(bin_edges)
# bin_centers = bin_edges[2:end]
inds = zeros(Int, length(df_iag.depth))
for i in eachindex(df_iag.depth)
    # find the index
    if df_iag.depth[i] < first(bin_edges)
        idx = firstindex(bin_edges)
    elseif df_iag.depth[i] > last(bin_edges)
        idx = lastindex(bin_edges)
    else
        idx = findlast(x -> x .<= df_iag.depth[i], bin_edges)
    end
    inds[i] = idx
end

# now bin the convective blueshifts
n_in_bin = zeros(length(bin_centers))
blueshift_avg = zeros(length(bin_centers))
blueshift_med = zeros(length(bin_centers))
blueshift_std = zeros(length(bin_centers))
blueshift_err = zeros(length(bin_centers))
blueshift_wts = zeros(length(bin_centers))
for i in eachindex(bin_centers)
    data = filter(!isnan, df_iag.blueshift[inds .== i])
    n_in_bin[i] = length(data)
    blueshift_avg[i] = mean(data)
    blueshift_med[i] = median(data)
    blueshift_std[i] = std(data)
    blueshift_err[i] = blueshift_std[i] / sqrt(length(data))
    blueshift_wts[i] = 1.0 / (sum((data .- blueshift_avg[i]).^2) / (n_in_bin[i] - 1))
end

# fit a polynomial to the binned data
@. model(x,p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3
p0 = [-504.891, -43.7963, -145.560, 884.308] # from IAG paper Reiners et al.
# blueshift_wts = 1.0 ./  blueshift_std.^2
blueshift_wts = ones(length(bin_centers))   # doesn't seem like Reiners et al. do WLS
pfit = curve_fit(model, bin_centers, blueshift_med, blueshift_wts, p0)
xs = range(bin_centers[1], bin_centers[end], length=1000)
ys = model(xs, pfit.param)

# estimate standard error on each param
if !all(blueshift_wts .== one(eltype(blueshift_wts)))
    se = stderror(pfit)

    # get uncertainty on fit to shade
    nsamps = 1000
    param_rand = zeros(nsamps, length(pfit.param))
    for i in eachindex(pfit.param)
        param_rand[:,i] = rand(Normal(pfit.param[i], se[i]), nsamps)
    end

    out = zeros(length(xs), nsamps)
    for i in 1:nsamps
        out[:, i] = model(xs, param_rand[i,:])
    end
    std_ys = dropdims(std(out, dims=2), dims=2)
    # plt.fill_between(xs, ys.-std_ys, y2=ys.+std_ys, color="k", alpha=0.25)
end

# get blueshifts to simulation
blueshifts = model(bin_centers, pfit.param)
rvs_avg, rvs_std = blueshift_vs_depth(bin_centers, blueshifts=[])
# rvs_avg, rvs_std = blueshift_vs_depth(bin_centers, blueshifts=zeros(length(bin_centers)))

# plot the result
fig, ax1 = plt.subplots()
ax1.scatter(depths, blueshift, s=1, alpha=0.5, c="k", label=L"{\rm IAG}")
ax1.errorbar(bin_centers, blueshift_med, yerr=blueshift_std, fmt=".", capsize=2.0,
             markersize=10, c="black", label=L"{\rm Binned\ IAG}")
plt.plot(xs, ys, "k--")
ax1.scatter(bin_centers, rvs_avg, c="tab:blue", s=25, label=L"{\rm Synthetic}", zorder=3)
plt.xlim(0.0, 1.0)
ax1.set_ylim(-1000, 400)
ax1.set_xlabel(L"{\rm Depth}")
ax1.set_ylabel(L"{\rm Convective\ Blueshift\ (ms^{-1})}")
ax1.legend(ncol=3, fontsize=11, loc="upper center")
fig.savefig(plotdir * "conv_blueshift.pdf")
plt.clf(); plt.close()
