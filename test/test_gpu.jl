using CUDA
using GRASS
using Test

@testset "GPU" begin
# values for SpecParams
lines = [5434.5]
depths = [0.25]
resolution = 7e5
spec = SpecParams(lines=lines, depths=depths, resolution=resolution)
Nλ = length(spec.lambdas)

lines = [5434.5]
depths = [0.5]
fname = GRASS.soldir * "FeI_5434.h5"
resolution = 7e5

# sort the input data for use on GPU
soldata = GRASS.SolarData(fname)
sorted_data = GRASS.sort_data_for_gpu(soldata)
disc_mu_cpu = sorted_data[1]
disc_ax_cpu = sorted_data[2]
lenall_cpu = sorted_data[3]
bisall_cpu = sorted_data[4]
intall_cpu = sorted_data[5]
widall_cpu = sorted_data[6]

# values for DiskParams
N = 132
Nt = 50
u1 = 0.4
u2 = 0.26
polex = 0.0
poley = 1.0
polez = 0.0
pole = (polex, poley, polez)
disk = DiskParams(N=N, Nt=Nt, pole=pole, u1=u1, u2=u2)
spec = SpecParams(lines=lines, depths=depths, templates=[fname], resolution=resolution)

@testset "Testing array sorting" begin
    # function to test that the len dimensions are correct
    function test_array_sizes(arr)
        for i in 1:size(arr, 3)
            idx = findfirst(row -> all(iszero.(row)), [arr[:,j,i] for j in 1:size(arr, 2)])
            if !isnothing(idx)
                @test lenall_cpu[i] == idx - 1
            else
                @test lenall_cpu[i] == size(arr, 2)
            end
        end
        return nothing
    end

    # now actually test len dims
    for x in sorted_data[4:end]
        test_array_sizes(x)
    end
end

@testset "Testing disk initialization" begin
    prec = Float64

    # move input data to gpu
    GRASS.@cusync begin
        disc_mu = CuArray{prec}(disc_mu_cpu)
        disc_ax = CuArray{Int32}(disc_ax_cpu)
        lenall_gpu = CuArray{Int32}(lenall_cpu)
        wavall_gpu = CuArray{prec}(wavall_cpu)
        widall_gpu = CuArray{prec}(widall_cpu)
        depall_gpu = CuArray{prec}(depall_cpu)
    end

    # allocate arrays for fresh copy of input data to copy to each loop
    GRASS.@cusync begin
        wavall_gpu_loop = CUDA.copy(wavall_gpu)
        depall_gpu_loop = CUDA.copy(depall_gpu)
    end

    # allocate memory for synthesis on the GPU
    GRASS.@cusync begin
        # indices, redshifts, and limb darkening
        tloop = CUDA.zeros(Int32, N, N)
        data_inds = CUDA.zeros(Int32, N, N)
        norm_terms = CUDA.zeros(prec, N, N)
        rot_shifts = CUDA.zeros(prec, N, N)

        # pre-allocated memory for interpolations
        starmap = CUDA.ones(prec, N, N, Nλ)
        lwavgrid = CUDA.zeros(prec, N, N, 100)
        rwavgrid = CUDA.zeros(prec, N, N, 100)
        allwavs = CUDA.zeros(prec, N, N, 200)
        allints = CUDA.zeros(prec, N, N, 200)
    end

    # move spatial and lambda grid to GPU
    GRASS.@cusync begin
        grid = CuArray{prec}(GRASS.make_grid(N))
        lambdas = CuArray{prec}(spec.lambdas)
    end
    grid_cpu = Array(grid)

    # set number of threads and blocks for N*N matrix gpu functions
    threads2 = (16, 16)
    blocks2 = cld(N^2, prod(threads2))

    # initialize values for data_inds, rot_shifts, and norm_terms
    GRASS.@cusync @cuda threads=threads2 blocks=blocks2 GRASS.initialize_arrays_for_gpu(data_inds, norm_terms, rot_shifts,
                                                                                  grid, disc_mu, disc_ax, u1,
                                                                                  u2, polex, poley, polez)

    # generate random indices for input data
    GRASS.@cusync @cuda threads=threads2 blocks=blocks2 GRASS.generate_tloop_gpu(tloop, grid, data_inds, lenall_gpu)

    # compute cpu values
    grid_cpu = GRASS.make_grid(N)
    norm_terms_cpu = zeros(N, N)
    vels_los_cpu = zeros(N, N)
    for i in eachindex(grid_cpu)
        for j in eachindex(grid_cpu)
            x = grid_cpu[i]
            y = grid_cpu[j]
            if GRASS.calc_r2(x,y) > 1.0
                continue
            end

            mu = GRASS.calc_mu(x,y)
            norm_terms_cpu[i,j] = GRASS.calc_norm_term(mu, N, u1, u2)
            vels_los_cpu[i,j] = GRASS.patch_velocity_los.(x, y, pole=pole)
        end
    end

        # compare cpu and gpu outputs
    @test all(isapprox.(norm_terms_cpu, Array(norm_terms), atol=1e-10))
    @test all(isapprox.(vels_los_cpu, Array(rot_shifts), atol=1e-10))


    # test that the correct input data is selected
    function get_gpu_data(idx, arr)
        idx2 = findfirst(row -> all(iszero.(row)), [arr[:,j,idx] for j in 1:size(arr, 2)])
        if isnothing(idx2)
            idx2 = size(arr,2) + 1
        end
        return arr[:,1:idx2-1,idx]
    end

    for i in eachindex(grid_cpu)
        for j in eachindex(grid_cpu)
            if GRASS.calc_r2(grid_cpu[i], grid_cpu[j]) > 1.0
                continue
            end

            # get cpu key
            key_cpu = GRASS.get_key_for_pos(grid_cpu[i], grid_cpu[j], unique(disc_mu_cpu), soldata.mu)

            # get gpu_key
            mu = GRASS.calc_mu(GRASS.calc_r2(grid_cpu[i], grid_cpu[j]))
            nn_mu_ind = GRASS.searchsortednearest_gpu(Array(disc_mu), mu)
            nn_ax_code = GRASS.find_nearest_ax_gpu(grid_cpu[i], grid_cpu[j])
            idx_gpu = GRASS.find_data_index_gpu(grid_cpu[i], grid_cpu[j], disc_mu, disc_ax)

            # now test the values
            @test soldata.wav[key_cpu] == get_gpu_data(idx_gpu, wavall_cpu)
            @test soldata.bis[key_cpu] == get_gpu_data(idx_gpu, bisall_cpu)
            @test soldata.dep[key_cpu] == get_gpu_data(idx_gpu, depall_cpu)
            @test soldata.wid[key_cpu] == get_gpu_data(idx_gpu, widall_cpu)
            @test Array(tloop)[i,j] <= lenall_cpu[idx_gpu]
        end
    end
end

@testset "Testing trimming" begin
    # set some trimming parameters
    depth = 0.25

    # move data to the GPU
    CUDA.@sync begin
        disc_mu = CuArray{Float64}(disc_mu_cpu)
        disc_ax = CuArray{Int32}(disc_ax_cpu)
        lenall_gpu = CuArray{Int32}(lenall_cpu)
        bisall_gpu = CuArray{Float64}(bisall_cpu)
        intall_gpu = CuArray{Float64}(intall_cpu)
        widall_gpu = CuArray{Float64}(widall_cpu)
    end

    # allocate memory for trim output on GPU
    CUDA.@sync begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
    end

    # loop over data for dik positions
    for n in eachindex(lenall_cpu)
        # loop over time slices for cpu
        for t in 1:lenall_cpu[n]
            # take time slice for disk position
            bist = view(bisall_cpu, :, t, n)
            intt = view(intall_cpu, :, t, n)
            widt = view(widall_cpu, :, t, n)

            # do the trim
            GRASS.trim_bisector!(depth, bist, intt)
        end

        # get slices on gpu
        CUDA.@sync begin
            bisall_gpu_in = CUDA.view(bisall_gpu, :, 1:lenall_cpu[n], n)
            intall_gpu_in = CUDA.view(intall_gpu, :, 1:lenall_cpu[n], n)

            # view of arrays to put modified bisectors in
            bisall_gpu_out = CUDA.view(bisall_gpu_loop, :, 1:lenall_cpu[n], n)
            intall_gpu_out = CUDA.view(intall_gpu_loop, :, 1:lenall_cpu[n], n)
        end

        # do the trim on the gpu
        threads1 = (16,16)
        blocks1 = cld(lenall_cpu[n] * 100, prod(threads1))
        CUDA.@sync @cuda threads=threads1 blocks=blocks1 GRASS.trim_bisector_gpu(depth, bisall_gpu_out,
                                                                           intall_gpu_out,
                                                                           bisall_gpu_in,
                                                                           intall_gpu_in)
    end

    for n in eachindex(lenall_cpu)
        plt.plot(Array(bisall_gpu)[:,:,n], Array(intall_gpu)[:,:,n], c="blue")
        plt.plot(Array(bisall_gpu_loop)[:,:,n], Array(intall_gpu_loop)[:,:,n], c="orange")
        plt.show()
    end

    # now test the trims match between cpu and gpu
    for n in eachindex(lenall_cpu)
        @test all(isapprox.(wavall_cpu[:,:,n], Array(wavall_gpu_loop)[:,:,n], atol=1e-10))
        @test all(isapprox.(bisall_cpu[:,:,n], Array(bisall_gpu_loop)[:,:,n], atol=1e-10))
        @test all(isapprox.(depall_cpu[:,:,n], Array(depall_gpu_loop)[:,:,n], atol=1e-10))
        @test all(isapprox.(widall_cpu[:,:,n], Array(widall_gpu_loop)[:,:,n], atol=1e-10))
    end
end

@testset "Testing synthesis" begin
    # set some trimming parameters
    top = NaN
    depth = 0.75

    # move data to the GPU
    CUDA.@sync begin
        lenall_gpu = CuArray(lenall_cpu)
        wavall_gpu = CuArray(wavall_cpu)
        bisall_gpu = CuArray(bisall_cpu)
        widall_gpu = CuArray(widall_cpu)
        depall_gpu = CuArray(depall_cpu)
    end

    # allocate memory for trim output on GPU
    CUDA.@sync begin
        wavall_gpu_loop = CUDA.copy(wavall_gpu)
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
        depall_gpu_loop = CUDA.copy(depall_gpu)
    end

    # loop over data for dik positions
    for n in eachindex(lenall_cpu)
        # loop over time slices for cpu
        for t in 1:lenall_cpu[n]
            # take time slice for disk position
            wavt = view(wavall_cpu, :, t, n)
            bist = view(bisall_cpu, :, t, n)
            dept = view(depall_cpu, :, t, n)
            widt = view(widall_cpu, :, t, n)

            # do the trim
            GRASS.trim_bisector_chop!(depth, wavt, bist, dept, widt, top=top)
        end

        # get slices on gpu
        CUDA.@sync begin
            wavt_gpu_out = CUDA.view(wavall_gpu_loop, :, 1:lenall_cpu[n], n)
            bist_gpu_out = CUDA.view(bisall_gpu_loop, :, 1:lenall_cpu[n], n)
            widt_gpu_out = CUDA.view(widall_gpu_loop, :, 1:lenall_cpu[n], n)
            dept_gpu_out = CUDA.view(depall_gpu_loop, :, 1:lenall_cpu[n], n)

            wavt_gpu_in = CUDA.view(wavall_gpu, :, 1:lenall_cpu[n], n)
            bist_gpu_in = CUDA.view(bisall_gpu, :, 1:lenall_cpu[n], n)
            widt_gpu_in = CUDA.view(widall_gpu, :, 1:lenall_cpu[n], n)
            dept_gpu_in = CUDA.view(depall_gpu, :, 1:lenall_cpu[n], n)
        end

        # do the trim on the gpu
        threads1 = (16,16)
        blocks1 = cld(lenall_cpu[n] * 100, prod(threads1))
        CUDA.@sync @cuda threads=threads1 blocks=blocks1 GRASS.trim_bisector_chop_gpu(depth,
                                                                                      wavt_gpu_out, bist_gpu_out,
                                                                                      dept_gpu_out, widt_gpu_out,
                                                                                      wavt_gpu_in, bist_gpu_in,
                                                                                      dept_gpu_in, widt_gpu_in, top)
    end

    disc_mu_cpu = disc_mu

    CUDA.@sync begin
        # move other data to the gpu
        grid = CuArray(GRASS.make_grid(N))
        disc_mu = CuArray(disc_mu)
        disc_ax = CuArray(disc_ax)
        lambdas = CuArray(spec.lambdas)

        # allocate memory for input data indices and normalization terms
        data_inds = CuArray{Int64,2}(undef, N, N)
        norm_terms = CuArray{Float64,2}(undef, N, N)
        rot_shifts = CuArray{Float64,2}(undef, N, N)
        λΔDs = CuArray{Float64,2}(undef, N, N)

        # allocate workspace memory on GPU
        starmap = CUDA.ones(Float64, N, N, Nλ)
        lwavgrid = CUDA.zeros(Float64, N, N, 100)
        rwavgrid = CUDA.zeros(Float64, N, N, 100)
        allwavs = CUDA.zeros(Float64, N, N, 200)
        allints = CUDA.zeros(Float64, N, N, 200)
    end


    # set number of threads and blocks for N*N matrix gpu functions
    threads2 = (16, 16)
    blocks2 = cld(N^2, prod(threads2))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads3 = (6,6,6)
    blocks3 = cld(N^2 * Nλ, prod(threads3))

    # get random starting indices and move it to GPU
    tloop = ones(Int64, N, N)
    tloop = CuArray(tloop)

    # initialize values for data_inds, tloop, and norm_terms
    CUDA.@sync @cuda threads=threads2 blocks=blocks2 GRASS.initialize_arrays_for_gpu(data_inds, norm_terms, rot_shifts,
                                                                                     grid, disc_mu, disc_ax, disk.u1,
                                                                                     disk.u2, polex, poley, polez)

    # fill workspace arrays on GPU
    threads4 = (6,6,6)
    blocks4 = cld(N^2 * 100, prod(threads4))
    CUDA.@sync @cuda threads=threads4 blocks=blocks4 GRASS.fill_workspace_arrays!(spec.lines[1], spec.conv_blueshifts[1],
                                                                                  grid, tloop, data_inds, rot_shifts,
                                                                                  λΔDs, wavall_gpu_loop,
                                                                                  bisall_gpu_loop, widall_gpu_loop,
                                                                                  depall_gpu_loop, lwavgrid,
                                                                                  rwavgrid, allwavs, allints)

    # fill workspace arrays on GPU
    threads4 = (6,6,6)
    blocks4 = cld(N^2 * 100, prod(threads4))
    CUDA.@sync @cuda threads=threads4 blocks=blocks4 GRASS.concatenate_workspace_arrays!(spec.lines[1], spec.depths[1],
                                                                                         spec.conv_blueshifts[1], grid, tloop,
                                                                                         data_inds, rot_shifts, λΔDs,
                                                                                         lenall_gpu, wavall_gpu_loop,
                                                                                         bisall_gpu_loop, widall_gpu_loop,
                                                                                         depall_gpu_loop, lwavgrid,
                                                                                         rwavgrid, allwavs, allints)

    # do the line synthesis
    CUDA.@sync @cuda threads=threads3 blocks=blocks3 GRASS.line_profile_gpu!(starmap, grid, lambdas, λΔDs, allwavs, allints)

    bad_ijs=[]

    # do stuff on cpu
    grid_cpu = GRASS.make_grid(N)
    for i in eachindex(grid_cpu)
        for j in eachindex(grid_cpu)
            # find position on disk and move to next iter if off disk
            x = grid_cpu[i]
            y = grid_cpu[j]
            r2 = GRASS.calc_r2(x, y)
            if r2 > 1.0
                continue
            end

            # find the nearest mu ind and ax code
            mu = GRASS.calc_mu(r2)
            nn_mu_ind = GRASS.searchsortednearest_gpu(disc_mu_cpu, mu)
            nn_ax_code = GRASS.find_nearest_ax_gpu(x, y)
            data_ind = GRASS.find_data_index_gpu(nn_mu_ind, nn_ax_code)

            lwavgrid_cpu = zeros(100)
            rwavgrid_cpu = zeros(100)
            allwavs_cpu = zeros(200)
            allints_cpu = zeros(200)
            prof_cpu = ones(length(spec.lambdas))

            # slice out bisectors for cpu
            wavt = view(wavall_cpu,:,1,data_ind)
            bist = view(bisall_cpu,:,1,data_ind)
            dept = view(depall_cpu,:,1,data_ind)
            widt = view(widall_cpu,:,1,data_ind)

            @test all(isapprox.(wavt, Array(wavall_gpu_loop)[:,1,data_ind]))
            @test all(isapprox.(bist, Array(bisall_gpu_loop)[:,1,data_ind]))
            @test all(isapprox.(dept, Array(depall_gpu_loop)[:,1,data_ind]))
            @test all(isapprox.(widt, Array(widall_gpu_loop)[:,1,data_ind]))


            # synthesize line on CPU
            GRASS.line_profile_cpu!(λΔDs[i,j], spec.lambdas, prof_cpu, wavt, dept, widt,
                                    lwavgrid_cpu, rwavgrid_cpu, allwavs_cpu, allints_cpu)

            # slice out the workspace arrays from gpu
            lwavgrid_gpu = Array(lwavgrid)[i,j,:]
            rwavgrid_gpu = Array(rwavgrid)[i,j,:]
            allwavs_gpu = Array(allwavs)[i,j,:]
            allints_gpu = Array(allints)[i,j,:]

            # now compare workspace values
            c1 = all(isapprox.(lwavgrid_cpu, Array(lwavgrid_gpu), atol=1e-10))
            c2 = all(isapprox.(rwavgrid_cpu, Array(rwavgrid_gpu), atol=1e-10))
            c3 = all(isapprox.(allwavs_cpu, Array(allwavs_gpu), atol=1e-10))
            c4 = all(isapprox.(allints_cpu, Array(allints_gpu), atol=1e-10))
            c5 = all(isapprox.(prof_cpu, Array(starmap)[i,j,:], atol=1e-10))

            if any([c1, c2, c3, c4, c5] .== false)
                push!(bad_ijs, (i,j))
                break
            end
        end
    end


end
"""

end
