using GRASS
using PyPlot
using LsqFit
using JLD2
using DataFrames
using CSV

Kostogryz_LD_file = DataFrame(CSV.File("data/Kostogryz_LD_300.csv", header = ["wavelength", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]))

function K_LD_best_fit(μ::T, wavelength::T) where T
    """
    limb darkening prescription for optical based on mu angle  
    """
    μ < zero(T) && return 0.0

    index = findmin(x->abs(x-wavelength), Kostogryz_LD_file[!, "wavelength"])[2]
    Kostogryz_LD_array = Kostogryz_LD_file[index, ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]]

    p0 = [1.0, 1.0]
    fit_K = curve_fit(GRASS.LD_model, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], collect(Kostogryz_LD_array), p0)
    p_opt_K = fit_K.param

    return GRASS.LD_model(μ, p_opt_K), p_opt_K[1], p_opt_K[2]
end

# mu array for figures
mu_arr = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
mu_zero_arr = range(0.0, step=0.05, stop=1.0)
# GRASS lines - A
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932, 15652.79]

best_u1 = Vector{Float64}(undef,size(wavelength)...)
best_u2 = Vector{Float64}(undef,size(wavelength)...)
best_u1_NL94 = Vector{Float64}(undef,size(wavelength)...)
best_u2_NL94 = Vector{Float64}(undef,size(wavelength)...)
#for each GRASS wavelength determine NL94 + Kostogryz LD at mu array
for lambda in 1:length(wavelength)
    NL_LD = Vector{Float64}(undef,size(mu_arr)...)
    NL_LD_full = Vector{Float64}(undef,size(mu_zero_arr)...)
    Kos_LD = Vector{Float64}(undef,size(mu_arr)...)
    Kos_LD_full = Vector{Float64}(undef,size(mu_zero_arr)...)

    # wavelength index for NL94
    index = findmin(x->abs(x-wavelength[lambda]/10.0), GRASS.lambda_nm)[2]

    #iterate through mu + collect corresponding NL94 + Kostogryz LD
    for i in 1:length(mu_arr)
        NL_LD[i] = GRASS.a0[index] + GRASS.a1[index]*mu_arr[i] + GRASS.a2[index]*mu_arr[i]^2 + GRASS.a3[index]*mu_arr[i]^3 + GRASS.a4[index]*mu_arr[i]^4 + GRASS.a5[index]*mu_arr[i]^5
        Kos_LD[i], best_u1[lambda], best_u2[lambda] = K_LD_best_fit(mu_arr[i], wavelength[lambda] ./ 10)
    end

    p0 = [1.0, 1.0]
    fit_NL = curve_fit(GRASS.LD_model, mu_arr, NL_LD, p0)
    p_opt_NL = fit_NL.param
    best_u1_NL94[lambda] = p_opt_NL[1]
    best_u2_NL94[lambda] = p_opt_NL[2]
  
    for i in 1:length(mu_zero_arr)
        Kos_LD_full[i] = GRASS.LD_model(mu_zero_arr[i], [best_u1[lambda], best_u2[lambda]])
        NL_LD_full[i] = GRASS.LD_model(mu_zero_arr[i], p_opt_NL)
    end

    #figures
    plt.scatter(mu_arr, NL_LD, label = "NL94", color = "b")
    plt.scatter(mu_arr, Kos_LD, label = "Kostogryz", color = "r")
    plt.plot(mu_zero_arr, NL_LD_full, color = "b")
    plt.plot(mu_zero_arr, Kos_LD_full, color = "r")
    plt.title(wavelength[lambda])
    plt.legend()
    plt.text(mu_arr[4], 0.5, "NL94 - GRASS wavelength $(round(GRASS.lambda_nm[index] - wavelength[lambda]/10.0; digits = 3)) nm")
    index = findmin(x->abs(x-wavelength[lambda] ./ 10), Kostogryz_LD_file[!, "wavelength"])[2]
    plt.text(mu_arr[4], 0.55, "Kostogryz - GRASS wavelength $(round(collect(Kostogryz_LD_file[!, "wavelength"])[index] - wavelength[lambda]/10.0; digits = 3)) nm")
    plt.xlabel("mu")
    plt.ylabel("relative intensity")
    plt.savefig("eclipse_figures/LD_fit/300/$(wavelength[lambda]).png")
    plt.clf()
end

df = DataFrame()
df[!, "wavelength"] = wavelength 
df[!, "u1"] = best_u1_NL94
df[!, "u2"] = best_u2_NL94
CSV.write("data/quad_ld_coeff_NL94_new.csv", df)