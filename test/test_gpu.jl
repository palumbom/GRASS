# CPU/GPU parity: disk_sim and disk_sim_gpu are independent, so each oracles the other.
# seed_rng=true routes tloop generation through the CPU for both, and parity assumes both
# geometry precomputes assign the same disk-position keys. The GPU is not bitwise
# reproducible -- line_profile_gpu! accumulates prof with atomics -- hence rtol, not ==.

const PARITY_RTOL = 1e-12
const PARITY_N = 50    # coarse; parity does not depend on geometry
const PARITY_NT = 4
const FLOAT32_RTOL = 1e-2   # elementwise; the single-precision error concentrates in the line core

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

# The trim buffers are reused across lines within one disk_sim_gpu call but reallocated
# per template, so lines-per-template and the template loop are separate axes. Cross them:
# testing each at the other's minimum leaves the interaction untested.
@testset "Testing parity with several lines per template, several templates" begin
    a, b = parity_λof("FeI_5250.2"), parity_λof("FeI_5250.6")
    same_group = SpecParams(lines=[a-0.15, a+0.15, b-0.15, b+0.15], depths=[0.5,0.9,0.6,0.4],
                            templates=["FeI_5250.2","FeI_5250.2","FeI_5250.6","FeI_5250.6"])
    diff_group = SpecParams(lines=[λ5434-0.15, λ5434+0.15,
                                   parity_λof("FeI_5576")-0.15, parity_λof("FeI_5576")+0.15],
                            depths=[0.5,0.9,0.6,0.4],
                            templates=["FeI_5434","FeI_5434","FeI_5576","FeI_5576"],
                            resolution=1e5)
    # crossed with variability as well, since the width buffer has the same lifetime
    mixed_var = SpecParams(lines=[a-0.15, a+0.15, b-0.15, b+0.15], depths=[0.5,0.5,0.5,0.5],
                           templates=["FeI_5250.2","FeI_5250.2","FeI_5250.6","FeI_5250.6"],
                           variability=[true,false,true,false])
    for spec in (same_group, diff_group, mixed_var)
        fc, fg = parity_fluxes(spec, disk)
        @test parity_relerr(fg, fc) <= PARITY_RTOL
        @test parity_rv_delta(spec, fc, fg) < 1e-6
    end
end

# three lines in one template: the line loop is where the trim runs, so line count is a
# distinct axis from the two-line case above
@testset "Testing parity with three lines in one template" begin
    spec = SpecParams(lines=[5434.0, 5434.5, 5435.0], depths=[0.4, 0.6, 0.9],
                      templates=["FeI_5434", "FeI_5434", "FeI_5434"])
    fc, fg = parity_fluxes(spec, disk)
    @test parity_relerr(fg, fc) <= PARITY_RTOL
    @test parity_rv_delta(spec, fc, fg) < 1e-6
end

# Every test above passes seed_rng=true, which routes tloop generation through the CPU.
# generate_tloop_gpu! and the unseeded branch of _synth_gpu are otherwise never executed.
# Unseeded output is not reproducible, so assert structure rather than values.
@testset "Testing unseeded GPU synthesis" begin
    spec = SpecParams(lines=[λ5434], depths=[0.75], templates=["FeI_5434"])
    _, f = synthesize_spectra(spec, disk, seed_rng=false, use_gpu=true,
                              verbose=false, show_progress=false)
    @test size(f) == (length(spec.lambdas), PARITY_NT)
    @test all(isfinite, f)
    @test all(isapprox.(maximum(f, dims=1), 1.0, atol=1e-8))
    @test all(0.0 .< minimum(f, dims=1) .< 1.0 - 0.5)

    # two lines in one template, unseeded: exercises the line loop on that path too
    spec2 = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.9],
                       templates=["FeI_5434", "FeI_5434"])
    _, f2 = synthesize_spectra(spec2, disk, seed_rng=false, use_gpu=true,
                               verbose=false, show_progress=false)
    @test all(isfinite, f2)
    @test all(isapprox.(maximum(f2, dims=1), 1.0, atol=1e-8))
end

# trim_bisector_gpu! does its scalar arithmetic in eltype(intt_in); a Float64 depth or
# literal would promote the kernel to double precision under precision=Float32
@testset "Testing the Float32 precision path" begin
    spec = SpecParams(lines=[λ5434], depths=[0.75], templates=["FeI_5434"])
    _, f64 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=true, precision=Float64,
                                verbose=false, show_progress=false)
    _, f32 = synthesize_spectra(spec, disk, seed_rng=true, use_gpu=true, precision=Float32,
                                verbose=false, show_progress=false)
    @test all(isfinite, f32)
    # single precision overshoots the continuum by a few Float32 eps
    @test all(isapprox.(maximum(f32, dims=1), 1.0, atol=1e-5))
    # elementwise, not norm-based: a norm test is loosest exactly where the error lives
    @test parity_relerr(f32, f64) <= FLOAT32_RTOL

    # depth 0.9 scales and 0.5 chops, so both int_top branches run in one call
    mixed = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.9],
                       templates=["FeI_5434", "FeI_5434"])
    _, m64 = synthesize_spectra(mixed, disk, seed_rng=true, use_gpu=true, precision=Float64,
                                verbose=false, show_progress=false)
    _, m32 = synthesize_spectra(mixed, disk, seed_rng=true, use_gpu=true, precision=Float32,
                                verbose=false, show_progress=false)
    @test all(isfinite, m32)
    @test parity_relerr(m32, m64) <= FLOAT32_RTOL

    # the trim kernel's typed IR for Float32 arrays must contain no Float64 and no Union
    sol32 = GRASS.GPUSolarData(SolarData(fname=joinpath(GRASS.soldir, "FeI_5434.h5")), precision=Float32)
    io = IOBuffer()
    CUDA.@device_code_warntype io=io @cuda launch=false GRASS.trim_bisector_gpu!(Float32(0.75), true,
        sol32.dep_contrast, sol32.len, copy(sol32.bis), copy(sol32.int), copy(sol32.wid),
        sol32.bis, sol32.int, sol32.wid)
    ir = String(take!(io))
    @test !occursin("Float64", ir)
    @test !occursin("Union{Float32, Float64}", ir)
end

end
