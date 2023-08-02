@testset "Geometry" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :make_grid)
    @test isdefined(GRASS, :calc_r2)
    @test isdefined(GRASS, :calc_mu)
    @test isdefined(GRASS, :find_nearest_ax)
end

@testset "Testing back/forth conversions" begin
    t = (0.75, 0.5)
    @test (0.75^2 + 0.5^2) == GRASS.calc_r2(t)
    @test sqrt(1.0 - GRASS.calc_r2(t)) == GRASS.calc_mu(t)
    @test GRASS.make_grid(132) == GRASS.make_grid(N=132)
    @test GRASS.calc_mu(t) == GRASS.calc_mu(t[1], t[2])
    @test GRASS.calc_mu((0.0, 0.0)) == 1.0
    @test GRASS.calc_mu((1.0, 2.5)) == 0.0
end

@testset "Testing off-disk coordinate" begin
    t = (1.25, 1.25)
    @test iszero(GRASS.calc_mu(t))
    @test GRASS.calc_r2(t) > one(eltype(t))
end

@testset "Testing normalization" begin
    N = 132
    u1 = 0.0
    u2 = 0.0
    @test GRASS.calc_norm_term(0.0, 0.0, N, u1, u2) == pi ./ (2.0 * N^2)

    # test summation
    norm = 0
    for x in GRASS.make_grid(N)
        for y in GRASS.make_grid(N)
            norm += GRASS.calc_norm_term(x, y, N, 0.4, 0.26)
        end
    end
    @test isapprox(norm, 1.0, atol=1e-3)
end

end
