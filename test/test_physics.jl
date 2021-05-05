using Test

@testset "Physics" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :absorption_line)
    @test isdefined(GRASS, :lorentzian_line)
    @test isdefined(GRASS, :width_thermal)
    @test isdefined(GRASS, :quad_limb_darkening)
    @test isdefined(GRASS, :calc_intensity_normalization)
    @test isdefined(GRASS, :rotation_period)
    @test isdefined(GRASS, :patch_velocity_los)
end

@testset "Testing limb darkening" begin
    # disk center
    t1 = (0.0, 0.0)
    u1 = 0.4
    u2 = 0.26
    mu = GRASS.calc_mu(t1)
    @test GRASS.quad_limb_darkening(mu, u1, u2) == 1.0

    # off disk
    t2 = (1.25, 1.25)
    mu = GRASS.calc_mu(t2)
    @test GRASS.quad_limb_darkening(mu, u1, u2) == 0.0

    # no LD
    u1 = 0.0
    u2 = 0.0
    mu = range(0.1, 1.0, length=10)
    @test all(GRASS.quad_limb_darkening.(mu, u1, u2) .== 1.0)
end

@testset "Testing rotational velocity" begin
    # default pole
    @test GRASS.patch_velocity_los(0.0, 0.0) == 0.0
    @test GRASS.patch_velocity_los(1.0, 0.0) < 0.0
    @test GRASS.patch_velocity_los(-1.0, 0.0) > 0.0

    # pole at face
    pole = (0.0, 0.0, 1.0)
    @test GRASS.patch_velocity_los(0.0, 0.0, pole=pole) == 0.0
    @test GRASS.patch_velocity_los(1.0, 0.0, pole=pole) == 0.0
    @test GRASS.patch_velocity_los(-1.0, 0.0, pole=pole) == 0.0
end

@testset "Testing Gaussian line" begin
    xs = range(5434.0, 5435.0, step=0.01)
    ys = GRASS.absorption_line.(xs, mid=5434.5, width=0.1, depth=1.0)
    @test minimum(ys) ≈ 0.0 atol=1e-5
    @test maximum(ys) ≈ 1.0 atol=1e-5
end

end
