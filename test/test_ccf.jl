using Test
using Random
using Statistics
import EchelleCCFs

@testset "CCF" begin

@testset "Testing function definitions" begin
    @test isdefined(SS, :calc_ccf)
    @test isdefined(SS, :calc_rvs_from_ccf)
    @test isdefined(EchelleCCFs, :calc_doppler_factor)
end

@testset "Testing velocity injection/extraction" begin
    # get line centers
    vel = 0.5 .* sin.(range(0, 2π, length=50)) # 50 m/s sin curve
    ΔλD = EchelleCCFs.calc_doppler_factor.(vel)
    λvs = ΔλD .* 5434.5

    # measure velocities for each Doppler shifted line
    res = 1e6
    xs = range(5433.0, 5436.0, step=5434.5/res)
    ccfs = []
    v_grid = []
    for i in eachindex(λvs)
        ys = SS.absorption_line.(xs, mid=λvs[i], width=0.1, depth=1.0)
        v_grid, ccf = calc_ccf(xs, ys, [5434.5], [1.0], res, normalize=true)
        push!(ccfs, ccf)
    end
    ccfs = cat(ccfs..., dims=2)
    rvs, sigs = SS.calc_rvs_from_ccf(v_grid, ccfs)
    delta = vel .- rvs

    # now actually test
    @test maximum(delta) < 0.1 # 0.1 m/s = 10 cm/s
end

end
