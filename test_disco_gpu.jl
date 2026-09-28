using GRASS
using CUDA

const DISCOJulia = GRASS.Eclipse.DISCOJulia
const DISCO_LINES = ("Fe5250", "Fe6152", "Fe6173")
const DISCO_DATA_DIR = get(
    ENV,
    "DISCO_PARAMS_DIR",
    normpath(joinpath(@__DIR__, "..", "DISCO", "GRASS", "data")),
)

function check(condition, description)
    condition || error("FAIL: $(description)")
    println("PASS: ", description)
end

function host_params_device(p::DISCOJulia.DISCOParams)
    DISCOJulia.DISCOParamsDevice(
        p.wavelength, p.mu_min, p.mu_max, Int32(p.n_mu),
        Int32(p.n_wave), Int32(p.n_pca),
        p.mean_profiles[1], p.mean_profiles[2], p.mean_profiles[3],
        p.eigenprofiles[1], p.eigenprofiles[2], p.eigenprofiles[3],
        p.pca_coeff_grids[1], p.pca_coeff_grids[2], p.pca_coeff_grids[3],
        p.e1_alpha, p.e1_zeta, p.e1_omega, p.e1_med,
        p.e2_alpha, p.e2_zeta, p.e2_omega, p.e2_med,
    )
end

function expected_disco_at_wavelength(params, mu, lambda_angstrom, z_rot, patch_id)
    wavelength = params.wavelength
    rest_nm = (lambda_angstrom / 10.0f0) / (1.0f0 + z_rot)
    if rest_nm < wavelength[1] || rest_nm > wavelength[end]
        return 1.0f0
    end

    index = min(searchsortedlast(wavelength, rest_nm), length(wavelength) - 1)
    fraction = (rest_nm - wavelength[index]) /
               (wavelength[index + 1] - wavelength[index])
    f0 = DISCOJulia.disco_intensity(
        params, mu, Int32(index), Int32(patch_id), UInt32(0))
    f1 = DISCOJulia.disco_intensity(
        params, mu, Int32(index + 1), Int32(patch_id), UInt32(0))
    return f0 + fraction * (f1 - f0)
end

function run_gpu_checks(line, model_path)
    println("Loading $(line) model from ", model_path)
    params = DISCOJulia.DISCOParams(model_path)
    params_gpu = DISCOJulia.to_device(params)
    params_host = host_params_device(params)

    mus = Float32[
        (params.mu_min + params.mu_max) / 2f0,
        params.mu_min + 0.8f0 * (params.mu_max - params.mu_min),
    ]
    dA_values = Float32[2.0, 1.0]
    ld_values = Float32[0.5, 1.0]
    ext_values = Float32[0.25, 0.75]
    z_rot_values = Float32[0.0, 1.0f-4]

    model_wavelength = params.wavelength
    sample_indices = unique(Int[
        1,
        2,
        cld(params.n_wave, 3),
        cld(2 * params.n_wave, 3),
        params.n_wave - 1,
        params.n_wave,
    ])
    wavelengths = Float32[model_wavelength[i] * 10f0 for i in sample_indices]
    pushfirst!(wavelengths, model_wavelength[1] * 10f0 - 0.02f0)
    push!(wavelengths, model_wavelength[end] * 10f0 + 0.02f0)

    for T in (Float32, Float64)
        mu_gpu = CUDA.CuArray{T}(reshape(mus, 2, 1))
        dA = CUDA.CuArray{T}(reshape(dA_values, 2, 1))
        ld = CUDA.CuArray{T}(reshape(ld_values, 2, 1, 1))
        ext = CUDA.CuArray{T}(reshape(ext_values, 2, 1, 1))
        z_rot = CUDA.CuArray{T}(reshape(z_rot_values, 2, 1))
        contrast = CUDA.zeros(T, 2, 1)
        lambda_gpu = CUDA.CuArray{T}(wavelengths)

        for ext_toggle in (0.0f0, 1.0f0)
            profile = CUDA.zeros(T, length(wavelengths))
            @cuda threads=(16, 16) blocks=(1, 1) GRASS.Eclipse.line_profile_disco_gpu!(
                1, profile, mu_gpu, ld, dA, ext, lambda_gpu, z_rot, contrast,
                ext_toggle, params_gpu, UInt32(0))

            weights = T[
                dA_values[i] * ld_values[i] *
                (ext_toggle == 1.0f0 ? ext_values[i] : 1.0f0)
                for i in eachindex(mus)
            ]
            sum_weights = sum(weights)
            expected = T[
                sum(
                    weights[i] * expected_disco_at_wavelength(
                        params_host, mus[i], lambda, z_rot_values[i], i)
                    for i in eachindex(mus)
                ) / sum_weights
                for lambda in wavelengths
            ]

            flux = CUDA.ones(T, length(wavelengths), 1)
            @cuda threads=256 blocks=1 GRASS.apply_line!(1, profile, flux, sum_weights)
            actual = Array(flux)[:, 1]
            label = "$(line), precision=$(T), extinction=$(ext_toggle == 1f0)"
            check(all(isfinite, actual), "finite flux ($label)")
            check(
                isapprox(actual, expected; rtol=2f-6, atol=2f-6),
                "GPU line synthesis matches CPU reference ($label); " *
                "max abs error=$(maximum(abs.(actual .- expected)))",
            )
        end
    end
end

CUDA.functional() || error(
    "CUDA is not functional on this node. Run this script in a GPU Slurm job.")

println("Running DISCO GPU line-kernel checks on ", CUDA.device())
for line in DISCO_LINES
    model_path = joinpath(DISCO_DATA_DIR, "disco_$(line)_params.h5")
    isfile(model_path) || error(
        "DISCO model file not found: $(model_path). Set DISCO_PARAMS_DIR to its directory.")
    run_gpu_checks(line, model_path)
end

println("All DISCO GPU checks passed for: ", join(DISCO_LINES, ", "))
