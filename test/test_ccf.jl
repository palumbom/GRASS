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

end
