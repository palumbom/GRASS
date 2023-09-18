@testset "Physics" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :absorption_line)
    @test isdefined(GRASS, :gaussian_line)
    @test isdefined(GRASS, :lorentzian_line)
    @test isdefined(GRASS, :width_thermal)
    @test isdefined(GRASS, :quad_limb_darkening)
    @test isdefined(GRASS, :rotation_period)
    @test isdefined(GRASS, :patch_velocity_los)
end

@testset "Testing limb darkening" begin
    # disk center
    u1 = 0.4
    u2 = 0.26
    @test GRASS.quad_limb_darkening(1.0, u1, u2) == 1.0

    # off disk
    @test GRASS.quad_limb_darkening(-0.25, u1, u2) == 0.0

    # no LD
    u1 = 0.0
    u2 = 0.0
    mu = range(0.1, 1.0, length=10)
    @test all(GRASS.quad_limb_darkening.(mu, u1, u2) .== 1.0)
end

@testset "Testing rotational velocity" begin
    # set up disk params
    disk90 = DiskParams(N=50, Nt=10, Nsubgrid=25, inclination=90.0)
    disk00 = DiskParams(N=50, Nt=10, Nsubgrid=25, inclination=0.0)

    # test equator on
    @test isapprox(GRASS.patch_velocity_los(deg2rad(0.0), deg2rad(0.0), disk90) * 3e8, 0.0, atol=1e-8)
    @test isapprox(GRASS.patch_velocity_los(deg2rad(0.0), deg2rad(-90.0), disk90) * 3e8, 2069.2, atol=1e-1)
    @test isapprox(GRASS.patch_velocity_los(deg2rad(0.0), deg2rad(90.0), disk90) * 3e8, -2069.2, atol=1e-1)

    # test pole on
    @test isapprox(GRASS.patch_velocity_los(deg2rad(90.0), deg2rad(0.0), disk00) * 3e8, 0.0, atol=1e-8)
    @test isapprox(GRASS.patch_velocity_los(deg2rad(0.0), deg2rad(180), disk00) * 3e8, 0.0, atol=1e-8)
    @test isapprox(GRASS.patch_velocity_los(deg2rad(0.0), deg2rad(0.0), disk00) * 3e8, 0.0, atol=1e-8)
end

@testset "Testing Gaussian line" begin
    xs = range(5434.0, 5435.0, step=0.01)
    ys = GRASS.gaussian_line.(xs, mid=5434.5, width=0.1, depth=1.0)
    @test minimum(ys) ≈ 0.0 atol=1e-5
    @test maximum(ys) ≈ 1.0 atol=1e-5

    ys = GRASS.gaussian_line.(xs, mid=5434.5, width=0.1, depth=0.5)
    @test minimum(ys) ≈ 0.5 atol=1e-5
    @test maximum(ys) ≈ 1.0 atol=1e-5
end

end
