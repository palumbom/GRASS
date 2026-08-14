# CPU/GPU parity for the disk-integrated synthesis path.
#
# disk_sim and disk_sim_gpu are independent implementations of one model, so each is an
# oracle for the other. With seed_rng=true -- which routes tloop generation through the
# CPU for both paths -- they agree to ~5e-15 in flux and 0.0000 m/s in RV when correct.
#
# TOLERANCE. The GPU is not bitwise reproducible run to run: line_profile_gpu!
# accumulates into prof with CUDA.@atomic, so the reduction order varies between
# launches. Measured spread over repeated runs is 4.2e-15 relative. PARITY_RTOL sits
# ~240x above that floor and ~1e10 below the smallest divergence these tests exist to
# catch (6.8e-4 in flux). Do not tighten it to bitwise equality -- the reproducibility
# testset below will fail first and tell you why. Do not widen it without re-measuring
# the floor on the hardware in question.
#
# Parity also assumes the CPU and GPU geometry precomputes assign identical disk-position
# keys (get_key_for_pos vs find_data_index_gpu). If these tests fail inexplicably after a
# change to either precompute, check the keys before suspecting the physics.

const PARITY_RTOL = 1e-12
const PARITY_N = 50    # coarse on purpose: parity is geometry-independent, and N=197 is slow
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

    # the CPU path is bitwise deterministic; the baseline harness relies on this
    _, c1 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=false,
                               verbose=false, show_progress=false)
    _, c2 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=false,
                               verbose=false, show_progress=false)
    @test c1 == c2

    # the GPU path is not, because prof is accumulated with atomics. It must still land
    # inside the parity tolerance -- if this fails, the tolerance is too tight for the
    # hardware, and every other test in this file is about to fail for the wrong reason.
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

# Two lines sharing one template put both lines inside a single disk_sim_gpu call, so
# they share bisall_gpu_loop. The GPU trim kernel used to skip writing the bisector on
# its scaling branch, leaving the previous line's chopped values in place (60 m/s).
# The scale branch fires on 0% of profiles at depth <= 0.6, 58% at 0.75, 100% at 0.9,
# so the mixed pair below is the trigger and the shallow pair is the control.
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

# widall_gpu_loop has the same cross-line lifetime as bisall_gpu_loop, and only the
# non-variability branch of the trim kernel writes it. A fixed-width line used to
# overwrite every epoch with epoch 1, contaminating any variable line that read it
# afterwards -- across lines and across time steps, so both orderings mattered (0.6 m/s).
@testset "Testing width reuse with mixed variability" begin
    for var in ([false, true], [true, false], [true, true], [false, false])
        spec = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.5],
                          templates=["FeI_5434", "FeI_5434"], variability=var)
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

# trim_bisector_chop! resamples up to maximum(intt); the GPU kernel hardcoded 1.0. These
# agree only for templates reaching the continuum. FeI_5382 tops out at 0.98542, and
# chops at depth 0.1 but scales at 0.4, giving a trigger and a control from one template.
@testset "Testing chop endpoint below the continuum" begin
    for dep in (0.1, 0.4)
        spec = SpecParams(lines=[λ5382], depths=[dep], templates=["FeI_5382"])
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

# disk_sim zeroes skipped epochs; disk_sim_gpu left them at the CUDA.ones initialization.
# simulate_observations divides binned flux by the count of unskipped epochs, so a
# leftover continuum silently inflates the result.
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

# Consecutive templates in the same line group reuse tloop_init; templates in different
# groups regenerate it. Both branches of that decision need covering.
@testset "Testing parity across multiple templates" begin
    same_group = SpecParams(lines=[parity_λof("FeI_5250.2"), parity_λof("FeI_5250.6")],
                            depths=[0.6, 0.4], templates=["FeI_5250.2", "FeI_5250.6"])
    # cross-group templates are >=140 A apart, so R=7e5 would make a ~10^5 point grid;
    # this case only exercises the regenerate branch, which is grid-spacing independent
    diff_group = SpecParams(lines=[λ5434, parity_λof("FeI_5576")], depths=[0.6, 0.4],
                            templates=["FeI_5434", "FeI_5576"], resolution=1e5)
    for spec in (same_group, diff_group)
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

end
