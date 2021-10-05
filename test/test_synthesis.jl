using DataFrames


@testset "Synthesis" begin

# get data for use in tests
data = GRASS.SolarData()
wavt = data.wav[(:c, :mu10)][:,1]
bist = data.bis[(:c, :mu10)][:,1]
dept = data.dep[(:c, :mu10)][:,1]
widt = data.wid[(:c, :mu10)][:,1]

# set trimming parameters
top = NaN
dep = 0.75

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :trim_bisector_chop!)
    @test isdefined(GRASS, :line_profile!)
    @test isdefined(GRASS, :synthesize_spectra)
end

@testset "Testing line trimming" begin
    # copy the data
    wavt1 = copy(wavt)
    bist1 = copy(bist)
    dept1 = copy(dept)
    widt1 = copy(widt)

    # do the trimming
    GRASS.trim_bisector_chop!(dep, wavt1, bist1, dept1, widt1, top=top)

    # test that it wasn't changed
    @test wavt1 !== wavt
    @test bist1 !== bist
    @test dept1 !== dept
    @test widt1 !== widt

    # top and bottom
    @test bist1[1] == 1.0 - dep
    @test bist1[end] == 1.0
    @test wavt1[1] == wavt[findfirst(bist .>= 1.0 - dep)]
end

@testset "Testing line synthesis" begin
    # copy the data
    wavt1 = copy(wavt)
    bist1 = copy(bist)
    dept1 = copy(dept)
    widt1 = copy(widt)

    # do the trimming
    GRASS.trim_bisector_chop!(dep, wavt1, bist1, dept1, widt1, top=top)

    # allocate memory for line
    λs = range(5434.0, 5435.0, length=1000)
    lwavgrid = zeros(100)
    rwavgrid = zeros(100)
    allwavs = zeros(200)
    allints = zeros(200)
    prof = ones(length(λs))

    # synthesize line
    mid = 5434.5
    GRASS.line_profile!(mid, λs, prof, wavt1, dept1, widt1, lwavgrid, rwavgrid, allwavs, allints)

    # test it
    @test !all(prof .== 1.0)
    @test prof[1] == 1.0
    @test prof[end] == 1.0
    @test minimum(prof) == 1.0 - dep
    @test maximum(prof) == 1.0
end

@testset "Testing bisector in/out" begin
    # read in data
    wavt1 = copy(wavt)
    bist1 = copy(bist)
    dept1 = copy(dept)
    widt1 = copy(widt)

    # do the trimming
    GRASS.trim_bisector_chop!(dep, wavt1, bist1, dept1, widt1, top=top)

    # allocate memory for line
    λs = range(5434.0, 5435.0, length=1000)
    lwavgrid = zeros(100)
    rwavgrid = zeros(100)
    allwavs = zeros(200)
    allints = zeros(200)
    prof = ones(length(λs))

    # synthesize line
    mid = 5434.5
    GRASS.line_profile!(mid, λs, prof, wavt1, dept1, widt1, lwavgrid, rwavgrid, allwavs, allints)

    # measure the bisector of synth line
    wavo, biso = GRASS.measure_bisector(λs, prof, interpolate=false)
    wavo .-= mid

    # test that output bisector matches input bisector
    @test all(isapprox.(wavo[2:99], wavt1[2:99], atol=1e-3))
end

@testset "Testing disk-integrated spectrum synthesis" begin
    # set params
    N = 132
    Nt = 10
    res = 7e5

    # simulate spectra
    spec = SpecParams(lines=[5434.5], depths=[dep], resolution=res)
    disk = DiskParams(N=N, Nt=Nt)
    wavs, flux = synthesize_spectra(spec, disk)

    # ensure spectra have correct properties
    @test size(flux,2) == Nt
    @test maximum(flux[:,1]) == 1.0
    # @test isapprox(minimum(flux[:,1]), 1.0 - dep, atol=1e-1) # TODO: fix?

    # get velocities
    v_grid, ccf1 = calc_ccf(wavs, flux, spec, normalize=true)
    rvs, sigs = calc_rvs_from_ccf(v_grid, ccf1)

    # test the velocities
    @test calc_rms(rvs) < 1.0
end

end
