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

# get plot directories
grassdir, plotdir, datadir = check_plot_dirs()

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

# evaluate the model on a fine grid of depths
depths_out = range(0.01, 1.0, step=0.01)
blueshifts_out = model(depths_out, pfit.param)

# do nearest neighbor interpolation on std
sigmas_out = zeros(length(depths_out))
for i in eachindex(depths_out)
    idx = GRASS.searchsortednearest(bin_centers, depths_out[i])
    sigmas_out[i] = blueshift_std[idx]
end

# write the IAG data to disk
df_out = DataFrame("depth" => depths_out, "blueshift" => blueshifts_out, "sigma" => sigmas_out)
CSV.write(GRASS.datdir * "convective_blueshift.dat", df_out)
