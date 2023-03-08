function disk_sim_rm_gpu(spec::SpecParams, disk::DiskParams, planet::Planet,
                         soldata::SolarData, outspec::AA{T,2}; precision::String="double",
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
    # @cusync begin
    @cusync begin
        disc_mu = CuArray{precision}(disc_mu_cpu)
        disc_ax = CuArray{Int}(disc_ax_cpu)
        lenall_gpu = CuArray{Int}(lenall_cpu)
        wavall_gpu = CuArray{precision}(wavall_cpu)
        bisall_gpu = CuArray{precision}(bisall_cpu)
        widall_gpu = CuArray{precision}(widall_cpu)
        depall_gpu = CuArray{precision}(depall_cpu)
    end

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        wavall_gpu_loop = CUDA.copy(wavall_gpu)
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
        depall_gpu_loop = CUDA.copy(depall_gpu)
    end

    # allocate memory for synthesis on the GPU
    @cusync begin
        starmap = CUDA.ones(precision, N, N, Nλ)
        lwavgrid = CUDA.zeros(precision, N, N, 100)
        rwavgrid = CUDA.zeros(precision, N, N, 100)
        allwavs = CUDA.zeros(precision, N, N, 200)
        allints = CUDA.zeros(precision, N, N, 200)
    end

    # move other data to the gpu
    grid_cpu = make_grid_range(N)
    @cusync begin
        grid = CuArray{precision}(grid_cpu)
        lambdas = CuArray{precision}(spec.lambdas)
    end

    # allocate memory for input data indices and normalization terms
    @cusync begin
        data_inds = CuArray{Int32,2}(undef, N, N)
        norm_terms = CuArray{precision,2}(undef, N, N)
        rot_shifts = CuArray{precision,2}(undef, N, N)
        λΔDs = CuArray{precision,2}(undef, N, N)
    end

    # get planet positions
    tvec = get_simulation_times(disk)
    xpos, ypos = calc_planet_position(tvec, planet)

    # set number of threads and blocks for N*N matrix gpu functions
    # threads2 = (16, 16)
    threads2 = (16, 16)
    blocks2 = cld(N^2, prod(threads2))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    # threads3 = (7,7,7)
    threads3 = (7,7,7)
    blocks3 = cld(N^2 * Nλ, prod(threads3))

    # set number of threads and blocks for N*N*100 matrix gpu functions
    threads4 = (6,6,6)
    # threads4 = (4,4,8)
    blocks4 = cld(N^2 * 100, prod(threads4))

    # initialize values for data_inds, tloop, and norm_terms
    @cusync @cuda threads=threads2 blocks=blocks2 initialize_arrays_for_gpu(data_inds, norm_terms, rot_shifts,
                                                                            grid, disc_mu, disc_ax, u1,
                                                                            u2, polex, poley, polez)

    # initialize values of tloop on CPU
    # TODO: on GPU generate float, multiply by what I want, and then floor it
    data_inds_cpu = Array(data_inds)
    tloop = zeros(Int32, N, N)
    generate_indices_for_gpu(tloop, grid_cpu, data_inds_cpu, lenall_cpu, verbose=verbose, seed_rng=seed_rng)

    # move tloop to GPU
    @cusync tloop = CuArray{Int32}(tloop)

    # loop over time
    for t in 1:Nt
        # don't do all this work if skip_times is true
        if skip_times[t]
            continue
        end

        # slice out position of planet
        xplanet = xpos[t]
        yplanet = ypos[t]

        # initialize starmap with fresh copy of weights, accounting for planet position
        @cusync @captured @cuda calc_norm_term_rm_gpu(starmap, norm_terms, rot_shifts, grid, xplanet,
                                                      yplanet, planet.radius, u1, u2, polex, poley, polez)
        # @cusync starmap .= norm_terms

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            # pre-trim the data, loop over all disk positions
            for n in eachindex(lenall_cpu)
                @cusync begin
                    # get correct index for position of input data
                    wavall_gpu_out = CUDA.view(wavall_gpu_loop, :, 1:lenall_cpu[n], n)
                    bisall_gpu_out = CUDA.view(bisall_gpu_loop, :, 1:lenall_cpu[n], n)
                    widall_gpu_out = CUDA.view(widall_gpu_loop, :, 1:lenall_cpu[n], n)
                    depall_gpu_out = CUDA.view(depall_gpu_loop, :, 1:lenall_cpu[n], n)

                    wavall_gpu_in = CUDA.view(wavall_gpu, :, 1:lenall_cpu[n], n) .* spec.variability[l]
                    bisall_gpu_in = CUDA.view(bisall_gpu, :, 1:lenall_cpu[n], n)
                    widall_gpu_in = CUDA.view(widall_gpu, :, 1:lenall_cpu[n], n)
                    depall_gpu_in = CUDA.view(depall_gpu, :, 1:lenall_cpu[n], n)
                end

                # do the trim
                threads1 = (16,16)
                blocks1 = cld(lenall_cpu[n] * 100, prod(threads1))
                @cusync @captured @cuda threads=threads1 blocks=blocks1 trim_bisector_chop_gpu(spec.depths[l],
                                                                                               wavall_gpu_out, bisall_gpu_out,
                                                                                               depall_gpu_out, widall_gpu_out,
                                                                                               wavall_gpu_in, bisall_gpu_in,
                                                                                               depall_gpu_in, widall_gpu_in,
                                                                                               NaN)
            end

            # fill workspace arrays
            @cusync @captured @cuda threads=threads4 blocks=blocks4 fill_workspace_arrays!(spec.lines[l], spec.conv_blueshifts[l],
                                                                                           grid, tloop, data_inds, rot_shifts, λΔDs,
                                                                                           wavall_gpu_loop, widall_gpu_loop,
                                                                                           lwavgrid, rwavgrid)
            @cusync @captured @cuda threads=threads4 blocks=blocks4 concatenate_workspace_arrays!(grid, tloop, data_inds, depall_gpu_loop,
                                                                                                  lwavgrid, rwavgrid, allwavs, allints)

            # do the line synthesis
            @cusync @captured @cuda threads=threads3 blocks=blocks3 line_profile_gpu!(starmap, grid, lambdas, λΔDs, allwavs, allints)

            # do array reduction and move data from GPU to CPU
            @cusync outspec[:,t] .*= Array(CUDA.view(CUDA.sum(starmap, dims=(1,2)), 1, 1, :))

            # iterate tloop
            if t < Nt
                @cusync @cuda threads=threads2 blocks=blocks2 iterate_tloop_gpu(tloop, data_inds, lenall_gpu, grid)
            end
        end
    end
    CUDA.synchronize()
    return spec.lambdas, outspec
end

# TODO: allocate the most I'll need
# TODO: add up whole sun where i,j are normal pixels
# TODO: then calculate profiles for occulted cells and subtract off
function calc_norm_term_rm_gpu(star_map, norm_terms, rot_shifts, grid, xpos, ypos, r_planet, u1, u2, polex, poley, polez)
    # get indices from GPU blocks + threads
    idx = threadIdx().x + blockDim().x * (blockIdx().x-1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y-1)
    sdy = blockDim().y * gridDim().y

    # parallelized loop over grid
    for i in idx:sdx:CUDA.length(grid)
        for j in idy:sdy:CUDA.length(grid)
            @inbounds star_map[i,j,:] .= norm_terms[i,j]
        end
    end
    return nothing
end
