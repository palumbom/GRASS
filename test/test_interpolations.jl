@testset "Interpolation" begin

@testset "Testing in-place interpolation" begin
    # make Gaussian line
    x1 = range(5434.0, 5435.0, step=0.01)
    y1 = GRASS.gaussian_line.(x1, mid=5434.5, width=0.1, depth=1.0)

    # new grid to re-sample on
    x2 = range(5433.5, 5435.5, step=0.05)

    # test that the interpolation runs
    itp1 = GRASS.linear_interp(x1, y1, bc=1.0)

    # test boundary-conditions
    y2 = itp1.(x2)
    @test y2[1] == 1.0
    @test y2[end] == 1.0

    # test flat boundary condition
    itp2 = GRASS.linear_interp(x1, y1, bc=NaN)
    y2 = itp2.(x2)
    @test y1[1] == y2[1]
    @test y1[end] == y2[end]
end


end
