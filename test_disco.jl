using GRASS
using HDF5

const DISCOJulia = GRASS.Eclipse.DISCOJulia

function check(condition, description)
    condition || error("FAIL: $(description)")
    println("PASS: ", description)
end

function sample_disco_params()
    DISCOJulia.DISCOParams(
        Float32[525.0, 525.1, 525.2],
        0.25f0, 1.0f0, 3, 3, 1,
        (
            Float32[0.2, 0.3, 0.4],
            Float32[0.4, 0.5, 0.6],
            Float32[0.6, 0.7, 0.8],
        ),
        (
            reshape(Float32[0.1, 0.2, 0.3], 1, 3),
            reshape(Float32[0.2, 0.1, 0.1], 1, 3),
            reshape(Float32[0.3, 0.2, 0.1], 1, 3),
        ),
        (
            reshape(Float32[0.1, 0.2, 0.3], 1, 3),
            reshape(Float32[0.0, 0.1, 0.2], 1, 3),
            reshape(Float32[0.1, 0.1, 0.1], 1, 3),
        ),
        zeros(Float32, 3), Float32[0.2, 0.4, 0.6],
        fill(0.1f0, 3), Float32[0.2, 0.4, 0.6],
        zeros(Float32, 3), Float32[0.5, 1.0, 1.5],
        fill(0.2f0, 3), Float32[0.5, 1.0, 1.5],
    )
end

function write_sample_disco_h5(path; mu_grid=Float32[0.25, 0.625, 1.0])
    p = sample_disco_params()
    h5open(path, "w") do f
        a = attrs(f)
        a["mu_min"] = p.mu_min
        a["mu_max"] = p.mu_max
        a["n_mu"] = p.n_mu
        a["n_wave"] = p.n_wave
        a["n_pca"] = p.n_pca
        f["wavelength"] = p.wavelength
        f["mu_grid"] = mu_grid

        pca = create_group(f, "pca")
        for (i, name) in enumerate(("GT", "OGR", "IgL"))
            group = create_group(pca, name)
            group["mean_profile"] = p.mean_profiles[i]
            group["eigenprofiles"] = p.eigenprofiles[i]
            group["pca_coeff_grid"] = p.pca_coeff_grids[i]
        end

        distributions = create_group(f, "distributions")
        for (name, alpha, zeta, omega, median) in (
            ("epsilon1", p.e1_alpha, p.e1_zeta, p.e1_omega, p.e1_med),
            ("epsilon2", p.e2_alpha, p.e2_zeta, p.e2_omega, p.e2_med),
        )
            group = create_group(distributions, name)
            group["alpha"] = alpha
            group["zeta"] = zeta
            group["omega"] = omega
            group["median"] = median
        end
    end
end

function sample_disco_device(p=sample_disco_params())
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

function expect_loader_error(path)
    try
        DISCOJulia.DISCOParams(path)
    catch err
        err isa ErrorException || rethrow()
        println("PASS: loader rejects a non-uniform μ grid")
        return
    end
    error("FAIL: loader accepted a non-uniform μ grid")
end

println("Checking DISCO HDF5 loading and interpolation...")
params = sample_disco_params()
device = sample_disco_device(params)

mktemp() do path, io
    close(io)
    write_sample_disco_h5(path)
    loaded = DISCOJulia.DISCOParams(path)
    check(loaded.wavelength == params.wavelength, "HDF5 wavelength round-trip")
    check(loaded.mean_profiles == params.mean_profiles, "HDF5 component means")
    check(loaded.eigenprofiles == params.eigenprofiles, "HDF5 eigenprofiles")
    check(loaded.pca_coeff_grids == params.pca_coeff_grids, "HDF5 PCA coefficient grids")
    check(loaded.e1_med == params.e1_med && loaded.e2_med == params.e2_med,
          "HDF5 filling-factor medians")

    expected_gt = 0.2f0 + 0.15f0 * 0.1f0
    actual_gt = DISCOJulia.disco_component_value(
        device, 0.4375f0, Int32(1), Int32(1))
    check(isapprox(actual_gt, expected_gt), "linear μ interpolation of a PCA component")

    gt = expected_gt
    ogr = 0.4f0 + 0.05f0 * 0.2f0
    igl = 0.6f0 + 0.1f0 * 0.3f0
    e1, e2 = 0.3f0, 0.75f0
    expected = gt * e1 +
               ogr * (e2 * (1f0 - e1) / (1f0 + e2)) +
               igl * ((1f0 - e1) / (1f0 + e2))
    actual = DISCOJulia.disco_intensity(
        device, 0.4375f0, Int32(1), Int32(1), UInt32(0))
    check(isapprox(actual, expected), "deterministic median filling-factor mixture")

    low = DISCOJulia.disco_intensity(
        device, 0.1f0, Int32(1), Int32(1), UInt32(0))
    low_edge = DISCOJulia.disco_intensity(
        device, params.mu_min, Int32(1), Int32(1), UInt32(0))
    high = DISCOJulia.disco_intensity(
        device, 1.1f0, Int32(1), Int32(1), UInt32(0))
    high_edge = DISCOJulia.disco_intensity(
        device, params.mu_max, Int32(1), Int32(1), UInt32(0))
    check(isapprox(low, low_edge) && isapprox(high, high_edge),
          "μ values outside the model domain clamp to the endpoints")
end

mktemp() do path, io
    close(io)
    write_sample_disco_h5(path; mu_grid=Float32[0.25, 0.6, 1.0])
    expect_loader_error(path)
end

const DISCO_LINES = ("Fe5250", "Fe6152", "Fe6173")
const DISCO_DATA_DIR = get(
    ENV,
    "DISCO_PARAMS_DIR",
    normpath(joinpath(@__DIR__, "..", "DISCO", "GRASS", "data")),
)

for line in DISCO_LINES
    model_path = joinpath(DISCO_DATA_DIR, "disco_$(line)_params.h5")
    isfile(model_path) || error(
        "DISCO model file not found: $(model_path). Set DISCO_PARAMS_DIR to its directory.")

    model = DISCOJulia.DISCOParams(model_path)
    check(model.n_wave > 2 && model.n_pca > 0, "$(line): load model")
    check(length(model.wavelength) == model.n_wave, "$(line): wavelength shape")
    check(all(diff(model.wavelength) .> 0f0), "$(line): increasing wavelength grid")
    check(all(isfinite, model.mean_profiles[1]), "$(line): finite component profiles")

    line_device = DISCOJulia.DISCOParamsDevice(
        model.wavelength, model.mu_min, model.mu_max, Int32(model.n_mu),
        Int32(model.n_wave), Int32(model.n_pca),
        model.mean_profiles[1], model.mean_profiles[2], model.mean_profiles[3],
        model.eigenprofiles[1], model.eigenprofiles[2], model.eigenprofiles[3],
        model.pca_coeff_grids[1], model.pca_coeff_grids[2], model.pca_coeff_grids[3],
        model.e1_alpha, model.e1_zeta, model.e1_omega, model.e1_med,
        model.e2_alpha, model.e2_zeta, model.e2_omega, model.e2_med,
    )
    center_profile = Float32[
        DISCOJulia.disco_intensity(
            line_device, 1.0f0, Int32(wave), Int32(1), UInt32(0))
        for wave in 1:model.n_wave
    ]
    check(all(isfinite, center_profile), "$(line): finite synthesized profile")
end

println("All DISCO CPU checks passed for: ", join(DISCO_LINES, ", "))
