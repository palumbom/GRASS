using Test

@testset "Geometry" begin

@testset "Testing function definitions" begin
    @test isdefined(SS, :make_grid)
    @test isdefined(SS, :calc_r2)
    @test isdefined(SS, :calc_mu)
    @test isdefined(SS, :find_nearest_mu)
    @test isdefined(SS, :find_nearest_ax)
end

@testset "Testing back/forth conversions" begin
    t = (0.75, 0.5)
    @test (0.75^2 + 0.5^2) == SS.calc_r2(t)
    @test sqrt(1.0 - SS.calc_r2(t)) == SS.calc_mu(t)
    @test SS.make_grid(256) == SS.make_grid(N=256)
    @test SS.calc_mu(t) == SS.calc_mu(t[1], t[2])
    @test SS.calc_mu((0.0, 0.0)) == 1.0
    @test SS.calc_mu((1.0, 2.5)) == 0.0
end

@testset "Testing off-disk coordinate" begin
    t = (1.25, 1.25)
    @test iszero(SS.calc_mu(t))
    @test SS.calc_r2(t) > one(eltype(t))
end

@testset "Testing normalization" begin
    N = 50
    u1 = 0.0
    u2 = 0.0
    @test SS.norm_term(0.0, 0.0, N, u1, u2) == pi ./ (2.0 * N^2)
end

end
