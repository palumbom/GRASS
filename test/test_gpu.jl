# CPU/GPU parity: disk_sim and disk_sim_gpu are independent, so each oracles the other.
# seed_rng=true routes tloop generation through the CPU for both, and parity assumes both
# geometry precomputes assign the same disk-position keys. The GPU is not bitwise
# reproducible -- line_profile_gpu! accumulates prof with atomics -- hence rtol, not ==.

const PARITY_RTOL = 1e-12
const PARITY_N = 50    # coarse; parity does not depend on geometry
const PARITY_NT = 4

parity_relerr(a, b) = maximum(abs.(a .- b) ./ max.(abs.(b), eps()))
parity_λof(name) = GRASS.get_template_wavelength(joinpath(GRASS.soldir, name * ".h5"))

# returns (flux_cpu, flux_gpu) for the same inputs and the same seeded tloop
function parity_fluxes(spec, disk; skip=falses(disk.Nt))
    _, fc = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=false,
                               skip_times=skip, verbose=false, show_progress=false)
    _, fg = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=true,
                               skip_times=skip, verbose=false, show_progress=false)
    return fc, fg
end

function parity_rv_delta(spec, fc, fg)
    vc, ccfc = calc_ccf(spec.lambdas, fc, spec, normalize=true)
    vg, ccfg = calc_ccf(spec.lambdas, fg, spec, normalize=true)
    rc, _ = calc_rvs_from_ccf(vc, ccfc)
    rg, _ = calc_rvs_from_ccf(vg, ccfg)
    fin = .!isnan.(rc) .& .!isnan.(rg)
    return any(fin) ? maximum(abs.(rc[fin] .- rg[fin])) : 0.0
end

@testset "CPU/GPU parity" begin

disk = DiskParams(N=PARITY_N, Nt=PARITY_NT)
λ5434 = parity_λof("FeI_5434")
λ5382 = parity_λof("FeI_5382")

@testset "Testing reproducibility of each path" begin
    spec = SpecParams(lines=[λ5434], depths=[0.75], templates=["FeI_5434"])

    # the cpu path is bitwise deterministic
    _, c1 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=false,
                               verbose=false, show_progress=false)
    _, c2 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=false,
                               verbose=false, show_progress=false)
    @test c1 == c2

    # the gpu path is not, but must stay inside the parity tolerance
    _, g1 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=true,
                               verbose=false, show_progress=false)
    _, g2 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=true,
                               verbose=false, show_progress=false)
    @test isapprox(g1, g2, rtol=PARITY_RTOL)
end

@testset "Testing parity for a single line" begin
    spec = SpecParams(lines=[λ5434], depths=[0.75], templates=["FeI_5434"])
    fc, fg = parity_fluxes(spec, disk)
    @test size(fc) == size(fg)
    @test parity_relerr(fg, fc) <= PARITY_RTOL
    @test parity_rv_delta(spec, fc, fg) < 1e-6
end

# two lines in one template share bisall_gpu_loop; depth 0.9 scales, 0.5 chops
@testset "Testing bisector reuse across the line loop" begin
    trigger = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.9],
                         templates=["FeI_5434", "FeI_5434"])
    control = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.5],
                         templates=["FeI_5434", "FeI_5434"])
    for (spec, name) in ((trigger, "mixed chop/scale"), (control, "both chop"))
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

# widall_gpu_loop is shared too, and only the fixed-width branch writes it
@testset "Testing width reuse with mixed variability" begin
    for var in ([false, true], [true, false], [true, true], [false, false])
        spec = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.5],
                          templates=["FeI_5434", "FeI_5434"], variability=var)
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

# FeI_5382 intensities stop short of the continuum; it chops at 0.1 and scales at 0.4
@testset "Testing chop endpoint below the continuum" begin
    for dep in (0.1, 0.4)
        spec = SpecParams(lines=[λ5382], depths=[dep], templates=["FeI_5382"])
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

# skipped epochs must be zero on both paths; binning divides by the unskipped count
@testset "Testing skip_times parity and semantics" begin
    spec = SpecParams(lines=[λ5434], depths=[0.75], templates=["FeI_5434"])
    skip = falses(PARITY_NT); skip[2] = true; skip[3] = true
    fc, fg = parity_fluxes(spec, disk, skip=skip)

    @test all(iszero, fc[:, skip])
    @test all(iszero, fg[:, skip])
    @test !any(iszero, fc[:, .!skip])
    @test !any(iszero, fg[:, .!skip])
    @test parity_relerr(fg, fc) <= PARITY_RTOL
end

# same-group templates reuse tloop_init, different groups regenerate it
@testset "Testing parity across multiple templates" begin
    same_group = SpecParams(lines=[parity_λof("FeI_5250.2"), parity_λof("FeI_5250.6")],
                            depths=[0.6, 0.4], templates=["FeI_5250.2", "FeI_5250.6"])
    # cross-group lines are far apart; drop the resolution to keep the grid small
    diff_group = SpecParams(lines=[λ5434, parity_λof("FeI_5576")], depths=[0.6, 0.4],
                            templates=["FeI_5434", "FeI_5576"], resolution=1e5)
    for spec in (same_group, diff_group)
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

end
