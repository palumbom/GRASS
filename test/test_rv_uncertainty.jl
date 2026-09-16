using Random
using Statistics
import EchelleCCFs

@testset "RV uncertainty" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :measure_rv_from_ccf_gaussian)
    @test isdefined(GRASS, :ccf_vcov)
    @test isdefined(GRASS, :ccf_stderror)
    @test isdefined(GRASS, :est_full_width)
    @test isdefined(GRASS, :find_idx_at_and_around_minimum)
    @test isdefined(GRASS, :gaussian_line_helper)
end

# a synthetic inverted-Gaussian CCF standing in for a real one
function synthetic_ccf(; v0::Float64=0.0, depth::Float64=0.5,
                         width::Float64=3000.0, n::Int=301)
    vels = collect(range(-15e3, 15e3, length=n))
    ccf = 1.0 .- depth .* exp.(-0.5 .* ((vels .- v0) ./ width).^2)
    return vels, ccf
end

@testset "Testing kwarg defaults resolve" begin
    # est_full_width and find_idx_at_and_around_minimum default their kwargs to
    # unexported EchelleCCFs constants; calling without kwargs is what catches a
    # missing import, since Julia only evaluates a default when it is taken
    vels, ccf = synthetic_ccf()

    fw = GRASS.est_full_width(vels, ccf)
    @test isfinite(fw)
    @test fw > 0.0

    amin, inds = GRASS.find_idx_at_and_around_minimum(vels, ccf)
    @test amin isa Integer
    @test first(inds) >= 1
    @test last(inds) <= length(vels)
end

@testset "Testing est_full_width when no width can be found" begin
    # NaN in the ccf poisons extrema, so every `ccf .<= target_val` comparison
    # is false and findfirst returns nothing. this is the only way to reach that
    # branch: target_val is always >= minccf, and minccf is an element of ccf,
    # so a finite ccf always has at least one pixel at or below the target
    vels, ccf = synthetic_ccf()
    ccf[50] = NaN
    @test isnan(GRASS.est_full_width(vels, ccf))

    # a flat ccf has zero depth but a well-defined width: the full span
    flat = ones(length(vels))
    @test GRASS.est_full_width(vels, flat) == vels[end] - vels[1]
end

@testset "Testing gaussian fit recovers an injected velocity" begin
    v_inject = 250.0
    vels, ccf = synthetic_ccf(v0=v_inject)
    ccf_var = fill(1e-6, length(ccf))

    mrv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=0.75)
    out = GRASS.measure_rv_from_ccf_gaussian(vels, ccf, ccf_var, mrv)

    @test isfinite(out.rv)
    @test isfinite(out.σ_rv)
    @test out.σ_rv > 0.0
    @test abs(out.rv - v_inject) < 1.0  # m/s
end

@testset "Testing an all-zero ccf returns NaN" begin
    vels = collect(range(-15e3, 15e3, length=101))
    ccf = zeros(length(vels))
    ccf_var = fill(1e-6, length(ccf))

    mrv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=0.75)
    out = GRASS.measure_rv_from_ccf_gaussian(vels, ccf, ccf_var, mrv)

    @test isnan(out.rv)
    @test isnan(out.σ_rv)
end

@testset "Testing the input variance is not mutated" begin
    # zero-variance pixels are substituted with eps() internally; that
    # substitution must not reach the caller's array
    vels, ccf = synthetic_ccf()
    ccf_var = fill(1e-6, length(ccf))
    ccf_var[1] = 0.0
    ccf_var[end] = 0.0
    ccf_var_before = copy(ccf_var)

    mrv = EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=0.75)
    GRASS.measure_rv_from_ccf_gaussian(vels, ccf, ccf_var, mrv)

    @test ccf_var == ccf_var_before
end

end
