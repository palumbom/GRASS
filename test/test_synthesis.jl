using Test
using CSV
using DataFrames

@testset "Synthesis" begin

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :trim_bisector_chop!)
    @test isdefined(GRASS, :line_profile!)
    @test isdefined(GRASS, :line_from_bis!)
end

@testset "Testing line trimming" begin
    # read in data
    df = CSV.read(GRASS.moddir * "test/synth_data.csv", DataFrame)
    wavt = copy(df.wavt)
    bist = copy(df.bist)
    dept = copy(df.dept)
    widt = copy(df.widt)

    # do the trimming
    depth = 0.75
    top = 0.9
    GRASS.trim_bisector_chop!(depth, wavt, bist, dept, widt, top=top)

    # test that it wasn't changed
    @test wavt !== df.wavt
    @test bist !== df.bist
    @test dept !== df.dept
    @test widt !== df.widt

    # top and bottom
    @test bist[1] == 1.0 - depth
    @test bist[end] == 1.0
    @test wavt[1] == df.wavt[findfirst(df.bist .>= 1.0 - depth)]
end

@testset "Testing line synthesis" begin
    # read in data
    df = CSV.read(GRASS.moddir * "test/synth_data.csv", DataFrame)
    wavt = copy(df.wavt)
    bist = copy(df.bist)
    dept = copy(df.dept)
    widt = copy(df.widt)

    # do the trimming
    depth = 0.75
    top = 0.9
    GRASS.trim_bisector_chop!(depth, wavt, bist, dept, widt, top=top)

    # allocate memory for line
    λs = range(5434.0, 5435.0, length=1000)
    lwavgrid = zeros(100)
    rwavgrid = zeros(100)
    allwavs = zeros(200)
    allints = zeros(200)
    prof = ones(length(λs))

    # synthesize line
    mid = 5434.5
    dep = 0.75
    GRASS.line_from_bis!(mid, λs, prof, wavt, dept, widt, lwavgrid, rwavgrid, allwavs, allints)

    # test it
    @test !all(prof .== 1.0)
    @test prof[1] == 1.0
    @test prof[end] == 1.0
    @test minimum(prof) == 1.0 - depth
    @test maximum(prof) == 1.0
end

@testset "Testing bisector in/out" begin
    # read in data
    df = CSV.read(GRASS.moddir * "test/synth_data.csv", DataFrame)
    wavt = copy(df.wavt)
    bist = copy(df.bist)
    dept = copy(df.dept)
    widt = copy(df.widt)

    # do the trimming
    depth = 0.75
    top = 0.9
    GRASS.trim_bisector_chop!(depth, wavt, bist, dept, widt, top=top)

    # allocate memory for line
    λs = range(5434.0, 5435.0, length=1000)
    lwavgrid = zeros(100)
    rwavgrid = zeros(100)
    allwavs = zeros(200)
    allints = zeros(200)
    prof = ones(length(λs))

    # synthesize line
    mid = 5434.5
    GRASS.line_from_bis!(mid, λs, prof, wavt, dept, widt, lwavgrid, rwavgrid, allwavs, allints)

    # measure the bisector of synth line
    wavo, biso = GRASS.measure_bisector(λs, prof, interpolate=false)
    wavo .-= mid

    # test that output bisector matches input bisector
    @test all(isapprox.(wavo[2:99], wavt[2:99], atol=1e-3))
end

end
