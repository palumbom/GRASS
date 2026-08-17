using DataFrames

@testset "Synthesis" begin

# get data for use in tests
data = GRASS.SolarData()
bist = view(data.bis[(:c, :mu10)], :, 1)
intt = view(data.int[(:c, :mu10)], :, 1)
widt = view(data.wid[(:c, :mu10)], :, 1)

# set trimming parameters
dep = 0.75

@testset "Testing function definitions" begin
    @test isdefined(GRASS, :trim_bisector!)
    @test isdefined(GRASS, :line_profile_cpu!)
    @test isdefined(GRASS, :synthesize_spectra)
end

@testset "Testing line trimming" begin
    # copy the data
    bist1 = copy(bist)
    intt1 = copy(intt)
    widt1 = copy(widt)
    bist2 = copy(bist)
    intt2 = copy(intt)
    widt2 = copy(widt)

    # do the trimming
    dep1 = 0.75
    dep2 = 0.95
    GRASS.trim_bisector!(dep1, bist1, intt1)
    GRASS.trim_bisector!(dep2, bist2, intt2)

    # test that it wasn't changed
    @test bist1 !== bist
    @test intt1 !== intt
    @test widt1 !== widt
    @test bist2 !== bist
    @test intt2 !== intt
    @test widt2 !== widt
    @test bist1 !== bist2
    @test intt1 !== intt2
    @test widt1 !== widt2

    # top and bottom
    @test intt1[1] == 1.0 - dep1
    @test intt2[1] == 1.0 - dep2
    @test intt1[end] == 1.0
    @test intt2[end] == 1.0
end

@testset "Testing line synthesis" begin
    # copy the data
    bist1 = copy(bist)
    intt1 = copy(intt)
    widt1 = copy(widt)

    # do the trimming
    GRASS.trim_bisector!(dep, bist1, intt1)

    # allocate memory for line
    λs = range(5432.0, 5437.0, step=5434.5/7e5)
    lwavgrid = zeros(100)
    rwavgrid = zeros(100)
    allwavs = zeros(200)
    allints = zeros(200)
    prof = zeros(length(λs))
    weight = 1.0

    # synthesize line
    mid = 5434.5
    GRASS.line_profile_cpu!(mid, weight, λs, prof, bist1, intt1, widt1,
                            lwavgrid, rwavgrid, allwavs, allints)

    # test it
    @test !all(prof .== 1.0)
    @test prof[1] == 1.0
    @test prof[end] == 1.0
    @test isapprox(minimum(prof), 1.0 - dep, atol=1e-3)
    @test isapprox(maximum(prof), 1.0, atol=1e-3)
end

@testset "Testing bisector in/out" begin
    # copy the data
    bist1 = copy(bist)
    intt1 = copy(intt)
    widt1 = copy(widt)

    # do the trimming
    GRASS.trim_bisector!(dep, bist1, intt1)
    
    # allocate memory for line
    λs = range(5432.0, 5437.0, step=5434.5/7e5)
    lwavgrid = zeros(100)
    rwavgrid = zeros(100)
    allwavs = zeros(200)
    allints = zeros(200)
    prof = zeros(length(λs))
    weight = 1.0

    # synthesize line
    mid = 5434.5
    GRASS.line_profile_cpu!(mid, weight, λs, prof, bist1, intt1, widt1,
                            lwavgrid, rwavgrid, allwavs, allints)

    # measure the bisector of synth line
    biso, into = GRASS.calc_bisector(λs, prof, nflux=100)
    biso .-= mid

    # test that output bisector matches input bisector
    @test all(isapprox.(biso[2:99], bist1[2:99], atol=1e-2))
end

@testset "Testing disk-integrated spectrum synthesis" begin
    # set params
    Nt = 50
    res = 7e5

    # simulate spectra
    spec = SpecParams(lines=[5434.5], depths=[dep], resolution=res, templates=["FeI_5434"])
    disk = DiskParams(Nt=Nt)
    wavs, flux = synthesize_spectra(spec, disk, verbose=false)

    # ensure spectra have correct properties
    @test size(flux,2) == Nt
    @test all(isapprox.(maximum(flux, dims=1), 1.0, atol=1e-8))
    @test all(isapprox.(minimum(flux, dims=1), 1.0 - dep[1], atol=1e-1))

    # get velocities
    v_grid, ccf1 = calc_ccf(wavs, flux, spec, normalize=true)
    rvs, sigs = calc_rvs_from_ccf(v_grid, ccf1)

    # test the velocities
    @test calc_rms(rvs) < 1.0
end

# cpu-only tests; the parity suite in test_gpu.jl always skips in CI
@testset "Testing template grouping" begin
    # the synthesis driver uses this in a boolean context, so it must return Bool
    @test GRASS.in_same_group("FeI_5250.2", "FeI_5250.6") isa Bool
    @test GRASS.in_same_group("FeI_5250.2", "FeI_5250.6")
    @test !GRASS.in_same_group("FeI_5250.2", "FeI_6301")

    # unlisted templates match nothing
    @test GRASS.in_same_group("SomeLine_1234", "FeI_5250.6") isa Bool
    @test !GRASS.in_same_group("SomeLine_1234", "FeI_5250.6")
    @test !GRASS.in_same_group("FeI_5250.2", "SomeLine_1234")
end

@testset "Testing synthesis argument validation" begin
    Nt = 4
    spec = SpecParams(lines=[5434.5], depths=[dep], templates=["FeI_5434"])
    disk = DiskParams(N=50, Nt=Nt)

    # a wrong length is a BoundsError or silently wrong flux columns
    @test_throws AssertionError synthesize_spectra(spec, disk, skip_times=falses(Nt - 1),
                                                   verbose=false, show_progress=false)
    @test_throws AssertionError synthesize_spectra(spec, disk, skip_times=falses(Nt + 1),
                                                   verbose=false, show_progress=false)
    @test_throws AssertionError synthesize_spectra(spec, disk, precision=Int,
                                                   verbose=false, show_progress=false)
end

@testset "Testing skip_times semantics" begin
    Nt = 4
    spec = SpecParams(lines=[5434.5], depths=[dep], templates=["FeI_5434"])
    disk = DiskParams(N=50, Nt=Nt)
    skip = falses(Nt); skip[2] = true; skip[3] = true

    wavs, flux = synthesize_spectra(spec, disk, skip_times=skip, verbose=false,
                                    show_progress=false)

    # skipped epochs are exactly zero, not continuum; binning relies on this
    @test size(flux, 2) == Nt
    @test all(iszero, flux[:, skip])
    @test !any(iszero, flux[:, .!skip])
    @test all(isapprox.(maximum(flux[:, .!skip], dims=1), 1.0, atol=1e-8))
end

# the guard precedes the use_gpu branch, so this runs without a GPU
@testset "Testing resolved synthesis argument validation" begin
    Nt = 4
    spec = SpecParams(lines=[5434.5], depths=[dep], templates=["FeI_5434"])
    disk = DiskParams(N=50, Nt=Nt)
    μ_bins = [0.25, 0.55, 0.85]

    # a wrong length is a BoundsError or silently wrong flux columns
    @test_throws AssertionError GRASS.synthesize_spectra_resolved(μ_bins, spec, disk,
                                    skip_times=falses(Nt - 1), verbose=false,
                                    show_progress=false)
    @test_throws AssertionError GRASS.synthesize_spectra_resolved(μ_bins, spec, disk,
                                    skip_times=falses(Nt + 1), verbose=false,
                                    show_progress=false)
end

@testset "Testing multiple lines from one template" begin
    Nt = 2
    # two lines sharing a template exercise the line loop inside a single disk_sim call
    spec = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.5],
                      templates=["FeI_5434", "FeI_5434"])
    disk = DiskParams(N=50, Nt=Nt)
    wavs, flux = synthesize_spectra(spec, disk, verbose=false, show_progress=false)

    @test size(flux, 2) == Nt
    @test all(isapprox.(maximum(flux, dims=1), 1.0, atol=1e-8))
    @test all(minimum(flux, dims=1) .< 0.75)

    # mixed variability must synthesize too
    spec_mixed = SpecParams(lines=[5434.2, 5434.8], depths=[0.5, 0.5],
                            templates=["FeI_5434", "FeI_5434"], variability=[false, true])
    wavs2, flux2 = synthesize_spectra(spec_mixed, disk, verbose=false, show_progress=false)
    @test all(isapprox.(maximum(flux2, dims=1), 1.0, atol=1e-8))
    @test all(minimum(flux2, dims=1) .< 0.75)
end

end
