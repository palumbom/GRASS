@testset "Geometry" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :make_grid)
    @test isdefined(GRASS, :calc_mu)
    @test isdefined(GRASS, :calc_dA)
    @test isdefined(GRASS, :sphere_to_cart)
    @test isdefined(GRASS, :find_nearest_ax_code)
    @test isdefined(GRASS, :get_grid_centers)
    @test isdefined(GRASS, :DiskParams)
    @test isdefined(GRASS, :SynthWorkspace)
    @test isdefined(GRASS, :precompute_quantities!)
end

@testset "Testing back/forth conversions" begin
    # set vectors
    O⃗ = [0.0, 1e6, 0.0]

    @test isnan(GRASS.calc_mu([0.0, 0.0, 0.0], O⃗))
    @test isone(GRASS.calc_mu([0.0, 1.0, 0.0], O⃗))

    @test isapprox(GRASS.sphere_to_cart(1.0, 0.0, 0.0), [1.0, 0.0, 0.0])
    @test isapprox(GRASS.sphere_to_cart(1.0, deg2rad(90), 0.0), [0.0, 0.0, 1.0])
    @test isapprox(GRASS.sphere_to_cart(1.0, 0.0, deg2rad(90)), [0.0, 1.0, 0.0])

    xyz = GRASS.sphere_to_cart(1.0, deg2rad(45.0), deg2rad(45.0))
    @test sqrt(1.0 - (xyz[1]^2.0 + xyz[3]^2.0)) == GRASS.calc_mu(xyz, O⃗)
end

@testset "Testing off-disk coordinate" begin
    O⃗ = [0.0, 1e6, 0.0]
    @test GRASS.calc_mu([0.0, -1.0, 0.0], O⃗) == -1.0
    @test iszero(GRASS.calc_mu([1.0, 0.0, 0.0], O⃗))
    @test iszero(GRASS.calc_mu([2.0, 0.0, 0.0], O⃗))
    @test iszero(GRASS.calc_mu([0.0, 0.0, 1.0], O⃗))
end

@testset "Testing normalization" begin
    N = 50
    Nsubgrid = 20
    Nt = 1
    u1 = 0.0
    u2 = 0.0
<<<<<<< HEAD
    @test GRASS.calc_norm_term(0.0, 0.0, N, u1, u2) == pi ./ (2.0 * N^2)

    # test summation of normalization terms
    norm2 = sum(GRASS.calc_norm_terms(N, u1, u2))
    @test isapprox(norm2, 1.0, atol=1e-6)
=======

    disk = DiskParams(N=N, Nt=Nt, Nsubgrid=Nsubgrid, u1=u1, u2=u2)
    wsp = GRASS.SynthWorkspace(disk)

    @test isapprox(sum(wsp.wts), π, atol=1e-5)
>>>>>>> main
end

end
