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

# determine whether to use GPU
const use_gpu = CUDA.functional()

function blueshift_vs_depth(depths::AbstractArray{Float64,1}, use_gpu::Bool=use_gpu)
    # set up stuff for lines
    lines = [5434.5]
    templates = ["FeI_5434"]
    resolution = 7e5

    # allocate memory and loop over depths
    rvs = zeros(length(depths))
    for (idx, depth) in enumerate(depths)
        @printf(">>> Doing depth %.2f \r", depth)
        spec = SpecParams(lines=lines, depths=[depth], templates=templates, resolution=resolution, extrapolate=true)

        # synthesize spectra
        disk = DiskParams(N=132, Nt=100)
        lambdas1, outspec1 = synthesize_spectra(spec, disk, use_gpu=use_gpu)

        # calculate the RVs
        v_grid, ccf1 = calc_ccf(lambdas1, outspec1, spec, normalize=true)
        rvs1, sigs1 = calc_rvs_from_ccf(v_grid, ccf1)
        rvs[idx] = mean(rvs1)
    end
    return rvs
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
blueshift_wts = 1.0 ./  blueshift_std.^2
pfit = curve_fit(model, bin_centers, blueshift_med, blueshift_wts, p0)
xs = range(bin_centers[1], bin_centers[end], length=1000)
ys = model(xs, pfit.param)

# estimate standard error on each param
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

# now do the plot
function plot_iag_blueshift()
    # do plotting
    plt.scatter(depths, blueshift, c="k", alpha=0.25, s=1.0)
    plt.errorbar(bin_centers, blueshift_med, yerr=blueshift_std, c="k", fmt=".", capsize=2.0)
    plt.plot(xs, ys, "k--", label="Fit")
    # plt.fill_between(xs, ys.-std_ys, y2=ys.+std_ys, color="k", alpha=0.25)
    plt.xlim(0.0, 1.0)
    plt.ylim(-1000, 400)
    plt.xlabel(L"{\rm Line\ Depth}")
    plt.ylabel(L"{\rm Convective\ Blueshift\ (ms^{-1})}")
    # plt.savefig("iag_blueshift.pdf")
    plt.show()
    plt.clf(); plt.close()
    return nothing
end

# plot_iag_blueshift()

# get blueshifts inherent to simulation
rvs = blueshift_vs_depth(bin_centers)


# # write the IAG data to disk
# df_iag = DataFrame("depth" => iag_depths, "blueshift" => iag_blueshifts, "sigma" => iag_sigmas)
# CSV.write(desktop * "convective_blueshift_iag.dat", df_iag)
# CSV.write(desktop * "feI5434_input/convective_blueshift.dat", df_iag)

# # iterate until converged
# conv_diff = zeros(1001)
# for i in 1:3
#     depths1 = range(0.05, 0.95, step=0.05)
#     rvs1 = blueshift_vs_depth(depths1)

#     # fit the relationship
#     pfit1 = Polynomials.fit(depths1, mean(rvs1, dims=1)'[:,1], 3)
#     sim_depths = range(0.0, 1.0, length=1001)
#     sim_blueshifts = pfit1.(sim_depths)

#     # get old blueshifts
#     df_old = CSV.read(desktop * "convective_blueshift_iag.dat",
#                       DataFrame, header=1, delim=",", type=Float64)

#     # get the difference
#     this_diff = df_old.blueshift .- sim_blueshifts
#     conv_diff .+= this_diff

#     # write the corrected blueshifts to disk
#     df_new = DataFrame("depth" => sim_depths, "blueshift" => df_old.blueshift .+ conv_diff, "sigma" => df_old.sigma)
#     CSV.write(desktop * "feI5434_input/convective_blueshift.dat", df_new)
# end


