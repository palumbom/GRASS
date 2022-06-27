import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()

function disk_sim_gpu(spec::SpecParams, disk::DiskParams, soldata::SolarData,
                      outspec::AA{T,2}; precision::String="double",
                      seed_rng::Bool=false, verbose::Bool=false,
                      skip_times::BitVector=BitVector(zeros(disk.Nt))) where T<:AbstractFloat
    # set single or double precision
    if precision == "single"
        precision = Float32
        println(">>> Single precision is broken right now!!")
    elseif precision  == "double"
        precision = Float64
    else
        precision = Float64
    end

    # seeding
    if seed_rng
        println(" >>> Seeding currently not implemented on GPU!")
    end

    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # get pole component vectors and limb darkening parameters
    polex = convert(precision, disk.pole[1])
    poley = convert(precision, disk.pole[2])
    polez = convert(precision, disk.pole[3])
    u1 = convert(precision, disk.u1)
    u2 = convert(precision, disk.u2)

    # sort the input data for use on GPU
    sorted_data = sort_data_for_gpu(soldata)
    disc_mu_cpu = sorted_data[1]
    disc_ax_cpu = sorted_data[2]
    lenall_cpu = sorted_data[3]
    wavall_cpu = sorted_data[4]
    bisall_cpu = sorted_data[5]
    widall_cpu = sorted_data[6]
    depall_cpu = sorted_data[7]

    # move input data to gpu
    @cusync begin
        disc_mu = CuArray{precision}(disc_mu_cpu)
        disc_ax = CuArray{Int}(disc_ax_cpu)
        lenall_gpu = CuArray{Int}(lenall_cpu)
        wavall_gpu = CuArray{precision}(wavall_cpu)
        widall_gpu = CuArray{precision}(widall_cpu)
        depall_gpu = CuArray{precision}(depall_cpu)
    end

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        wavall_gpu_loop = CUDA.copy(wavall_gpu)
        depall_gpu_loop = CUDA.copy(depall_gpu)
    end

    # allocate memory for synthesis on the GPU
    outspec_temp = similar(outspec)
    @cusync begin
        # indices, redshifts, and limb darkening
        tloop = CUDA.zeros(Int32, N, N)
        data_inds = CUDA.zeros(Int32, N, N)
        norm_terms = CUDA.zeros(precision, N, N)
        rot_shifts = CUDA.zeros(precision, N, N)
        λΔDs = CUDA.zeros(precision, N, N)

        # pre-allocated memory for interpolations
        starmap = CUDA.ones(precision, N, N, Nλ)
        lwavgrid = CUDA.zeros(precision, N, N, 100)
        rwavgrid = CUDA.zeros(precision, N, N, 100)
        allwavs = CUDA.zeros(precision, N, N, 200)
        allints = CUDA.zeros(precision, N, N, 200)
    end

    # move spatial and lambda grid to GPU
    @cusync begin
        grid = CuArray{precision}(make_grid(N))
        lambdas = CuArray{precision}(spec.lambdas)
    end

    # set number of threads and blocks for trimming functions
    threads1 = (16,16)
    blocks1 = cld.(lenall_cpu .* 100, prod(threads1))

    # set number of threads and blocks for N*N matrix gpu functions
    threads2 = (16, 16)
    blocks2 = cld(N^2, prod(threads2))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads3 = (7,7,7)
    blocks3 = cld(N^2 * Nλ, prod(threads3))

    # set number of threads and blocks for N*N*100 matrix gpu functions
    threads4 = (6,6,6)
    blocks4 = cld(N^2 * 100, prod(threads4))

    # initialize values for data_inds, rot_shifts, and norm_terms
    @cusync @cuda threads=threads2 blocks=blocks2 initialize_arrays_for_gpu(data_inds, norm_terms, rot_shifts,
                                                                            grid, disc_mu, disc_ax, u1,
                                                                            u2, polex, poley, polez)

    # generate random indices for input data
    @cusync @cuda threads=threads2 blocks=blocks2 generate_tloop_gpu(tloop, grid, data_inds, lenall_gpu)

    # loop over time
    for t in 1:Nt
        # don't do all this work if skip_times is true
        if skip_times[t]
            continue
        end

        # initialize starmap with fresh copy of weights
        @cusync starmap .= norm_terms

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # pre-trim the data, loop over all disk positions
            for n in eachindex(lenall_cpu)
                epoch_range = 1:lenall_cpu[n]
                @cusync begin
                    # view of unaltered input data
                    wavall_gpu_in = CUDA.view(wavall_gpu, :, epoch_range, n) .* spec.variability[l]
                    depall_gpu_in = CUDA.view(depall_gpu, :, epoch_range, n)

                    # view of arrays to put modified bisectors in
                    wavall_gpu_out = CUDA.view(wavall_gpu_loop, :, epoch_range, n)
                    depall_gpu_out = CUDA.view(depall_gpu_loop, :, epoch_range, n)
                end

                # do the trim
                @cusync @captured @cuda threads=threads1 blocks=blocks1[n] trim_bisector_gpu(spec.depths[l],
                                                                                             wavall_gpu_out,
                                                                                             depall_gpu_out,
                                                                                             wavall_gpu_in,
                                                                                             depall_gpu_in)
            end

            # fill workspace arrays
            @cusync @captured @cuda threads=threads4 blocks=blocks4 fill_workspace_arrays!(spec.lines[l], spec.conv_blueshifts[l],
                                                                                           grid, tloop, data_inds, rot_shifts, λΔDs,
                                                                                           wavall_gpu_loop, widall_gpu,
                                                                                           lwavgrid, rwavgrid)
            @cusync @captured @cuda threads=threads4 blocks=blocks4 concatenate_workspace_arrays!(grid, tloop, data_inds, depall_gpu_loop,
                                                                                                  lwavgrid, rwavgrid, allwavs, allints)

            # do the line synthesis
            @cusync @captured @cuda threads=threads3 blocks=blocks3 line_profile_gpu!(starmap, grid, lambdas, λΔDs, allwavs, allints)

            # do array reduction and move data from GPU to CPU
            @cusync outspec[:,t] .*= Array(CUDA.view(CUDA.sum(starmap, dims=(1,2)), 1, 1, :))

            # iterate tloop
            @cusync @cuda threads=threads2 blocks=blocks2 iterate_tloop_gpu(tloop, data_inds, lenall_gpu, grid)
        end
    end
    # CUDA.synchronize()
    return spec.lambdas, outspec
end

# TODO: compare kernel launch cost to compute cost
# TODO: how many registers are needed?
