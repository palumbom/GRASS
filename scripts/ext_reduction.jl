using JLD2
using SPICE
using GRASS
using Revise
using CSV
using DataFrames
using Statistics
import PyPlot
plt = PyPlot
mpl = plt.matplotlib

function extinction(disk, neid_ext_coeff, idx1, idx3, zenith_mean, mu_grid, ext)
    for i in eachindex(disk.ϕc)
        for j in 1:disk.Nθ[i]
            if all(x -> x < 0.0, vec(mu_grid[i, j]))
                continue
            end

            boolean_mask = idx3[i,j] .== 1
            idx1_sum = sum(idx1[i,j])
            
            for ext_ind in 1:length(neid_ext_coeff)
                extin = map(x -> exp(-((1/cosd(x))*neid_ext_coeff[ext_ind])), zenith_mean[i,j])
                ext[i,j,ext_ind] = mean(view(extin, boolean_mask)) 
            end
        end
    end
end

function extinction_reduction(time, index)

    variable_names = ["dA_total_proj_mean", "zenith", "idx1", "idx3", "N", "mu_grid"]

    # Open the JLD2 file and read the variables into a dictionary
    geometry = jldopen("/storage/home/efg5335/work/GRASS/data/solar_disk/neid_october_exp_meter_N_50_$(index)000.jld2", "r") do file
        Dict(var => read(file, var) for var in variable_names)
    end
    
    dA_total_proj_mean = deepcopy(geometry["dA_total_proj_mean"])
    zenith_mean = deepcopy(geometry["zenith"])
    idx1 = deepcopy(geometry["idx1"])
    idx3 = deepcopy(geometry["idx3"])
    mu_grid = deepcopy(geometry["mu_grid"])

    intensity = jldopen("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Eclipse_Figures/ExposureMeter/data/neid_october_exp_meter_N_50_$(index)000_KSSD.jld2", "r") do file
        Dict(var => read(file, var) for var in ["mean_intensity"])
    end

    mean_intensity = deepcopy(intensity["mean_intensity"])

    # set up paramaters for disk
    N = geometry["N"]
    Nt = length(time)
    disk = GRASS.DiskParamsEclipse(N=N, Nt=Nt, Nsubgrid=10)
    ext_coeff = LinRange(0.09, 0.2, 20)
    wsp = GRASS.SynthWorkspaceEclipse(disk, length(ext_coeff), Nt, verbose=true)

    # get lines to construct templates
    lp = GRASS.LineProperties()
    name = GRASS.get_name(lp)
    λrest = GRASS.get_rest_wavelength(lp)
    lfile = GRASS.get_file(lp)

    intensity_ext_final = Vector{Vector{Vector{Float64}}}(undef,size(lp.λrest)...)

    # loop over lines
    Threads.@threads for i in eachindex(lp.λrest) 
        println("\t>>> Template: " * string(splitdir(lfile[i])[2]))

        intensity_ext_list = Vector{Vector{Float64}}(undef,size(time)...)
        for t in 1:Nt
            extinction(disk, ext_coeff, idx1[t], idx3[t], zenith_mean[t], mu_grid[t], wsp.ext)

            idx_grid = mean_intensity[i][t] .> 0.0

            inner_list = Vector{Float64}(undef,size(ext_coeff)...)
            for ext_ind in 1:length(ext_coeff)
                brightness = mean_intensity[i][t] .* dA_total_proj_mean[:, :, t] .* wsp.ext[:, :, ext_ind]
                inner_list[ext_ind] = sum(view(brightness, idx_grid))
            end
            intensity_ext_list[t] = inner_list
        end
        intensity_ext_final[i] = intensity_ext_list
    end

    return intensity_ext_final
end

#exposure meter data
exp_meter_data = DataFrame(CSV.File("data/exposure_meter_data.csv"))
exp_meter_time = exp_meter_data[!, "Array2"]

start = 1
end_value = 1000

intensity_ext_eclipse = Vector{Vector{Vector{Vector{Float64}}}}(undef,size(1:7)...)
Threads.@threads for j in 1:7
    if j == 7
        intensity_ext_eclipse[j] = extinction_reduction(exp_meter_time[start:length(exp_meter_time)], j)
    else
        intensity_ext_eclipse[j] = extinction_reduction(exp_meter_time[start:end_value], j)
    end

    global start += 1000
    global end_value = parse(Int, "$(j+1)000")
end

@save "extinction_chi_near1.jld2"
    jldopen("extinction_chi_near1.jld2", "a+") do file
        file["extinction_chi_data"] = deepcopy(intensity_ext_eclipse) 
end