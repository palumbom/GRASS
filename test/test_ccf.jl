using Random
using Statistics
import EchelleCCFs

@testset "CCF" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :fit_voigt)
    @test isdefined(GRASS, :calc_ccf)
    @test isdefined(GRASS, :calc_rvs_from_ccf)
    @test isdefined(EchelleCCFs, :calc_doppler_factor)
end

@testset "Testing velocity injection/extraction" begin
    # get line centers
    mid = 5434.5
    dep = 0.75
    amp = 0.5
    vel = amp .* sin.(range(0, 2π, length=100)) # in m/s
    ΔλD = EchelleCCFs.calc_doppler_factor.(vel)
    λvs = ΔλD .* mid

    # generate spectra
    res = 7e5
    buff = 1.25
    Δlnλ = (1.0 / res)
    wavs = exp.(range(log(mid-buff), log(mid+buff), step=Δlnλ))
    flux = zeros(length(wavs), length(λvs))
    for i in eachindex(λvs)
        flux[:,i] = GRASS.fit_voigt(wavs, [0.3, λvs[i], 0.04, 0.045])
    end

    # calculate ccf and measure velocity
    v_grid, ccf = calc_ccf(wavs, flux, [mid], [dep], res, normalize=true)
    rvs, sigs = GRASS.calc_rvs_from_ccf(v_grid, ccf)

    # subtract off mean
    rvs .-= mean(rvs)

    # test that the residuals are within 10 cm/s
    delta = vel .- rvs
    @test maximum(delta) < 0.01 # 0.01 m/s = 1 cm/s
end

@testset "Testing variance-weighted CCF and RVs" begin
    # single-epoch synthetic line, same construction as the injection test above
    mid = 5434.5
    dep = 0.75
    res = 7e5
    buff = 1.25
    Δlnλ = (1.0 / res)
    wavs = exp.(range(log(mid-buff), log(mid+buff), step=Δlnλ))
    flux = GRASS.fit_voigt(wavs, [0.3, mid, 0.04, 0.045])

    # constant per-pixel variance for a given continuum SNR
    snr = 500.0
    var = fill((1.0/snr)^2, length(wavs))

    # the variance path and the unnormalized plain path must agree exactly:
    # this is code motion, not a different algorithm
    v_grid, ccf, ccf_var = calc_ccf(wavs, flux, var, [mid], [dep], res)
    v_grid_plain, ccf_plain = calc_ccf(wavs, flux, [mid], [dep], res, normalize=false)
    @test v_grid == v_grid_plain
    @test ccf == ccf_plain

    @test length(ccf_var) == length(v_grid)
    @test all(ccf_var .>= 0.0)
    @test all(isfinite.(ccf_var))

    # the SpecParams wrapper must accept a Float64 variance, not just Float32.
    # checked with hasmethod rather than by constructing a SpecParams: its
    # constructor reads convective_blueshift.dat from datdir, and this file's
    # other testsets are deliberately free of any data dependency
    @test hasmethod(calc_ccf, (Vector{Float64}, Vector{Float64},
                               Vector{Float64}, GRASS.SpecParams{Float64}))
    @test hasmethod(calc_ccf, (Vector{Float64}, Vector{Float64},
                               Vector{Float32}, GRASS.SpecParams{Float64}))

    # the RV from the variance path lands on the same line core as the plain path
    out = calc_rvs_from_ccf(v_grid, ccf, ccf_var)
    rv_plain, _ = GRASS.calc_rvs_from_ccf(v_grid, ccf)
    @test isfinite(out.rv)
    @test isfinite(out.σ_rv)
    @test out.σ_rv > 0.0
    @test abs(out.rv - rv_plain) < 1.0  # m/s

    # the caller's variance array survives the call unchanged. the zero must be
    # planted: the internal substitution only writes where an element is exactly
    # 0.0, and a ccf variance off a tophat mask is strictly positive everywhere,
    # so without this the assertion holds against a mutating implementation too
    ccf_var[1] = 0.0
    ccf_var_before = copy(ccf_var)
    calc_rvs_from_ccf(v_grid, ccf, ccf_var)
    @test ccf_var == ccf_var_before
    @test ccf_var[1] == 0.0

    # normalize was silently ignored on the variance path; it is now rejected
    @test_throws MethodError calc_ccf(wavs, flux, var, [mid], [dep], res, normalize=true)

    # only the Gaussian path is implemented: a quadratic fit_type is rejected at the
    # boundary rather than failing on the missing init_guess_ccf_σ field mid-fit
    @test_throws AssertionError calc_rvs_from_ccf(v_grid, ccf, ccf_var,
                                                 fit_type=GRASS.QuadraticFit)
    @test calc_rvs_from_ccf(v_grid, ccf, ccf_var, fit_type=GRASS.GaussianFit).rv == out.rv
end

end
