@testset "Input Data" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :SpecParams)
    @test isdefined(GRASS, :DiskParams)
    @test isdefined(GRASS, :SolarData)
    @test isdefined(GRASS, :LineProperties)
end

@testset "Testing data read-in and sizes" begin
    # default read-in for 5434.5 line
    soldat = GRASS.SolarData()
    @test !isempty(soldat.bis)
    @test !isempty(soldat.int)
    @test !isempty(soldat.wid)
    @test keys(soldat.bis) == keys(soldat.int) == keys(soldat.wid)

    # test that len keyword is working
    k = [k for k in keys(soldat.len)]
    @test all(map(x -> soldat.len[x] == size(soldat.bis[x],2), k))
    @test all(map(x -> soldat.len[x] == size(soldat.int[x],2), k))
    @test all(map(x -> soldat.len[x] == size(soldat.wid[x],2), k))
end

end
