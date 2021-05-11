@testset "Input Data" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :SpecParams)
    @test isdefined(GRASS, :DiskParams)
    @test isdefined(GRASS, :SolarData)
end

@testset "Testing data read-in" begin
    soldat = GRASS.SolarData()
    @test !isempty(soldat.wav)
    @test !isempty(soldat.bis)
    @test length(keys(soldat.wav)) == 41
    @test keys(soldat.wav) == keys(soldat.bis)
    for k in keys(soldat.len)
        @test soldat.len[k] == size(soldat.wav[k], 2)
        @test soldat.len[k] == size(soldat.bis[k], 2)
        @test size(soldat.bis[k], 2) == size(soldat.wav[k], 2)
    end
end
