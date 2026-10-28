using PyCall
using CSV, DataFrames
using PyPlot
using HDF5
using CUDA
using Printf
using PyCall

py"""
import sys
sys.path.append('.')
"""

DISCO_ROOT = "/storage/home/efg5335/work/sw/DISCO"
include(joinpath(DISCO_ROOT, "line_creation/make_lines.jl"))
using .DISCOJulia

line = "Fe6173"
ref_wl = 6173.3344

function kernel_single_profile!(prof_d, p_device, mu, patch_id, epoch_seed)
    w = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if w <= length(prof_d)
        prof_d[w] = DISCOJulia.disco_intensity(
            p_device, mu, Int32(w), patch_id, epoch_seed
        )
    end
    return nothing
end

function get_disco_profile_gpu(
    p_device::DISCOJulia.DISCOParamsDevice,
    mu::Float32;
    patch_id::Int32=Int32(1),
    epoch_seed::UInt32=UInt32(0),
)
    n_wave = Int(p_device.n_wave)
    prof_d = CUDA.zeros(Float32, n_wave)
    
    threads = 256
    blocks = cld(n_wave, threads)
    
    @cuda threads=threads blocks=blocks kernel_single_profile!(prof_d, p_device, mu, patch_id, epoch_seed)
    CUDA.synchronize()
    
    return Array(prof_d) 
end

function average_bisector_velocity(profile, wavelengths, ref_wl; max_level=0.75, n_levels=100, c=2.99792458e8)
    I_cont = Float64(maximum(profile))
    core_idx = argmin(profile)
    I_min = Float64(profile[core_idx])
    I_top = Float64(max_level * I_cont)

    # Ensure valid range between core and max_level
    if I_top <= I_min
        I_top = I_min + 0.01 * (I_cont - I_min)
    end

    # Sample intensity levels uniformly from near the core up to I_top
    # (avoid exact minimum point to ensure well-defined left and right wings)
    intensities = range(I_min + 1e-5 * (I_top - I_min), I_top, length=n_levels)
    
    λ_bis_sum = 0.0
    n_valid = 0

    n_pts = length(profile)

    for I_k in intensities
        # 1. Find blue (left) wing crossing (index < core_idx)
        blue_idx = 0
        for i in core_idx:-1:2
            if profile[i-1] >= I_k >= profile[i]
                blue_idx = i - 1
                break
            end
        end

        # 2. Find red (right) wing crossing (index >= core_idx)
        red_idx = 0
        for j in core_idx:(n_pts - 1)
            if profile[j] <= I_k <= profile[j+1]
                red_idx = j
                break
            end
        end

        # Skip level if crossing is not bounded on both sides
        if blue_idx == 0 || red_idx == 0
            continue
        end

        # Linear interpolation on blue wing
        y1_b, y2_b = Float64(profile[blue_idx]), Float64(profile[blue_idx + 1])
        x1_b, x2_b = Float64(wavelengths[blue_idx]), Float64(wavelengths[blue_idx + 1])
        λ_blue = x1_b + (I_k - y1_b) * (x2_b - x1_b) / (y2_b - y1_b)

        # Linear interpolation on red wing
        y1_r, y2_r = Float64(profile[red_idx]), Float64(profile[red_idx + 1])
        x1_r, x2_r = Float64(wavelengths[red_idx]), Float64(wavelengths[red_idx + 1])
        λ_red = x1_r + (I_k - y1_r) * (x2_r - x1_r) / (y2_r - y1_r)

        # Bisector wavelength at intensity level I_k
        λ_bis_sum += 0.5 * (λ_blue + λ_red)
        n_valid += 1
    end

    if n_valid == 0
        # Fallback to discrete minimum if no levels could be sampled
        return c * (wavelengths[core_idx] - ref_wl) / ref_wl
    end

    # Average bisector wavelength from core up to 75% level
    mean_λ_bis = λ_bis_sum / n_valid

    # Convert average bisector wavelength shift to Doppler velocity in m/s
    return c * (mean_λ_bis - ref_wl) / ref_wl
end

function extract_cb(h5_path::String)
    p_host = DISCOJulia.DISCOParams(h5_path)
    p_device = DISCOJulia.to_device(p_host)

    MU_MIN = 0.0
    MU_MAX = 1.00
    N_MU_DENSE = 4096
    mus = range(MU_MIN, MU_MAX, N_MU_DENSE)

    wl_nm = p_host.wavelength
    wl_angstroms = wl_nm .* 10.0f0  

    jl_cb = Vector{Float64}(undef, length(mus))
    polynomial_cb = Vector{Float64}(undef, length(mus))
    original_cb = Union{Missing,Float64}[missing for _ in eachindex(mus)]

    reference_nm = ref_wl / 10
    for i in eachindex(mus)
        mu = Float32(mus[i])
        jl_prof = get_disco_profile_gpu(p_device, mu; epoch_seed=UInt32(0))
        jl_prof = jl_prof ./ maximum(jl_prof)

        jl_cb[i] = average_bisector_velocity(jl_prof, wl_angstroms, ref_wl; max_level=0.05, n_levels=100)
        polynomial_cb[i] =
            DISCOJulia.cb_polynomial(p_host, mu) / reference_nm * 3e8
        if mu >= p_host.mu_min
            original_cb[i] =
                DISCOJulia.original_cb_offset(p_host, mu) / reference_nm * 3e8
        end
    end

    df = DataFrame(
        Mu = Float64.(mus),
        DISCOCB = original_cb,
        PolynomialCB = polynomial_cb,
        MeasuredCB = jl_cb,
        constant = fill(p_host.cb_poly_c0_nm / reference_nm * 3e8, length(mus)),
        a = fill(
            p_host.cb_poly_c2_nm_per_mu2 / reference_nm * 3e8, length(mus)
        ),
        b = fill(
            p_host.cb_poly_c1_nm_per_mu / reference_nm * 3e8, length(mus)
        )
    )

    CSV.write("data/disco_$(line)_cb.csv", df)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    H5_PATH = "data/disco_$(line)_params.h5"
    extract_cb(H5_PATH)
end
