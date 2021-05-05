using Interpolations

@testset "Interpolation" begin

@testset "Testing in-place interpolation" begin
    # make Gaussian line
    x1 = range(5434.0, 5435.0, step=0.01)
    y1 = SS.absorption_line.(x1, mid=5434.5, width=0.1, depth=1.0)

    # new grid to re-sample on
    x2 = range(5433.5, 5435.5, step=0.05)

    # test that the interpolation runs
    itp1 = extrapolate(interpolate!(eltype(x1), (x1,), y1, Gridded(Linear())), 1.0)
    @test typeof(itp1) <: Interpolations.FilledExtrapolation

    # test boundary-conditions
    y2 = itp1.(x2)
    @test y2[1] == 1.0
    @test y2[end] == 1.0
end


end
