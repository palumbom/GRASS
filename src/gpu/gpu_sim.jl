function disk_sim_gpu(spec::SpecParams{T}, disk::DiskParams{T}, soldata::SolarData{T},
                      gpu_allocs::GPUAllocs, outspec::AA{T,2}; verbose::Bool=false,
                      seed_rng::Bool=false, precision::DataType=Float64,
                      skip_times::BitVector=falses(disk.Nt)) where T<:AF
    # get dimensions for memory alloc
    N = disk.N
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # convert scalars from disk params to desired precision
    ρs = convert(precision, disk.ρs)
    A = convert(precision, disk.A)
    B = convert(precision, disk.B)
    C = convert(precision, disk.C)
    u1 = convert(precision, disk.u1)
    u2 = convert(precision, disk.u2)

    # parse out composite type
    ϕc = gpu_allocs.ϕc
    θc = gpu_allocs.θc
    R_θ = gpu_allocs.R_θ
    O⃗ = gpu_allocs.O⃗
    λs = gpu_allocs.λs
    μs = gpu_allocs.μs
    tloop = gpu_allocs.tloop
    data_inds = gpu_allocs.data_inds
    norm_terms = gpu_allocs.norm_terms
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs
    starmap = gpu_allocs.starmap
    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints

    # sort the input data for use on GPU
    sorted_data = sort_data_for_gpu(soldata)
    disc_mu_cpu = sorted_data[1]
    disc_ax_cpu = sorted_data[2]
    lenall_cpu = sorted_data[3]
    cbsall_cpu = sorted_data[4]
    depcontrast_cpu = sorted_data[5]
    bisall_cpu = sorted_data[6]
    intall_cpu = sorted_data[7]
    widall_cpu = sorted_data[8]

    # set number of threads and blocks for N*N matrix gpu functions
    threads1 = (16,16)
    blocks1 = cld(N^2, prod(threads1))

    # set number of threads and blocks for trimming functions
    threads2 = (4,4,16)
    blocks2 = cld(length(lenall_cpu) * maximum(lenall_cpu) * 100, prod(threads2))

    # set number of threads and blocks for N*N*100 matrix gpu functions
    threads3 = (4,4,16)
    blocks3 = cld(N^2 * 100, prod(threads3))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads4 = (3,3,42)
    blocks4 = cld(N^2 * Nλ, prod(threads4))

    # move input data to gpu
    @cusync begin
        disc_mu_gpu = CuArray{precision}(disc_mu_cpu)
        disc_ax_gpu = CuArray{Int32}(disc_ax_cpu)
        lenall_gpu = CuArray{Int32}(lenall_cpu)
        cbsall_gpu = CuArray{precision}(cbsall_cpu)
        bisall_gpu = CuArray{precision}(bisall_cpu)
        intall_gpu = CuArray{precision}(intall_cpu)
        widall_gpu = CuArray{precision}(widall_cpu)
        depcontrast_gpu = CuArray{precision}(depcontrast_cpu)
    end

    # allocate arrays for fresh copy of input data to copy to each loop
    @cusync begin
        bisall_gpu_loop = CUDA.zeros(precision, CUDA.size(bisall_gpu))
        intall_gpu_loop = CUDA.zeros(precision, CUDA.size(intall_gpu))
        widall_gpu_loop = CUDA.zeros(precision, CUDA.size(widall_gpu))
    end

    # # # get launch parameters
    # # kernel = @cuda launch=false initialize_arrays_for_gpu(tloop, data_inds, norm_terms, z_rot, z_cbs,
    # #                                                       grid, disc_mu_gpu, disc_ax_gpu, lenall_gpu,
    # #                                                       cbsall_gpu, u1, u2, polex, poley, polez)

    # # config = launch_configuration(kernel.fun)
    # # threads1 = min(N^2, config.threads)
    # # blocks1 = cld(N^2, threads1)

    # # @show threads1
    # # @show blocks1

    # # kernel(tloop, data_inds, norm_terms, z_rot, z_cbs, grid,
    # #        disc_mu_gpu, disc_ax_gpu, lenall_gpu, cbsall_gpu,
    # #        u1, u2, polex, poley, polez; threads=threads1, blocks=blocks1)


    # initialize values for data_inds, tloop, dop_shifts, and norm_terms
    vec1 = CUDA.zeros(precision, length(ϕc), length(ϕc), 3)
    vec2 = CUDA.zeros(precision, length(ϕc), length(ϕc), 3)
    vec3 = CUDA.zeros(precision, length(ϕc), length(ϕc), 3)
    @cusync @cuda threads=threads1 blocks=blocks1 initialize_arrays_for_gpu(vec1, vec2, vec3, ρs, ϕc, θc, R_θ, O⃗, μs,
                                                                            tloop,
                                                                            data_inds, norm_terms,
                                                                            z_rot, z_cbs, disc_mu_gpu,
                                                                            disc_ax_gpu, lenall_gpu,
                                                                            cbsall_gpu, A, B, C, u1, u2)

    dat = Array(z_rot)
    # dat[dat .<= 0.0] .= NaN
    plt.pcolormesh(Array(vec1)[:,:,1], Array(vec1)[:,:,3], dat .* 3e8)
    plt.colorbar()

    # plopt = Array(z_rot)
    # cmap = mpl.cm.ScalarMappable(cmap="viridis")
    # clrs = cmap.to_rgba(plopt)
    # plt.plot_surface(Array(vec1[:,:,1]), Array(vec1[:,:,2]), Array(vec1[:,:,3]), facecolors=clrs)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.zlabel("z")
    plt.show()

    # # get weighted disk average cbs
    # @cusync sum_norm_terms = CUDA.sum(norm_terms)
    # @cusync z_cbs_avg = CUDA.sum(z_cbs .* norm_terms) / sum_norm_terms

    # # calculate how much extra shift is needed
    # extra_z = spec.conv_blueshifts .- z_cbs_avg

    # # loop over time
    # for t in 1:Nt
    #     # don't synthesize spectrum if skip_times is true, but iterate t index
    #     if skip_times[t]
    #         @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, data_inds, lenall_gpu, grid)
    #         continue
    #     end

    #     # initialize starmap with fresh copy of weights
    #     @cusync starmap .= norm_terms

    #     # loop over lines to synthesize
    #     for l in eachindex(spec.lines)
    #         @cusync begin
    #             CUDA.copyto!(bisall_gpu_loop, bisall_gpu)
    #             CUDA.copyto!(intall_gpu_loop, intall_gpu)
    #             CUDA.copyto!(widall_gpu_loop, widall_gpu)
    #         end

    #         # trim all the bisector data
    #         @cusync @captured @cuda threads=threads2 blocks=blocks2 trim_bisector_gpu!(spec.depths[l], spec.variability[l],
    #                                                                                    depcontrast_gpu, lenall_gpu,
    #                                                                                    bisall_gpu_loop, intall_gpu_loop,
    #                                                                                    widall_gpu_loop, bisall_gpu,
    #                                                                                    intall_gpu, widall_gpu)

    #         # assemble line shape on even int grid
    #         @cusync @captured @cuda threads=threads3 blocks=blocks3 fill_workspaces!(spec.lines[l], spec.variability[l],
    #                                                                                  extra_z[l], grid,
    #                                                                                  tloop, data_inds, z_rot, z_cbs,
    #                                                                                  bisall_gpu_loop, intall_gpu_loop,
    #                                                                                  widall_gpu_loop, allwavs, allints)

    #         # do the line synthesis, interp back onto wavelength grid
    #         @cusync @captured @cuda threads=threads4 blocks=blocks4 line_profile_gpu!(starmap, grid, λs, allwavs, allints)
    #     end

    #     # do array reduction and move data from GPU to CPU
    #     @cusync @inbounds outspec[:,t] .*= dropdims(Array(CUDA.sum(starmap, dims=(1,2))), dims=(1,2))

    #     # iterate tloop
    #     @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, data_inds, lenall_gpu, grid)
    # end

    # # ensure normalization
    # outspec ./= sum_norm_terms
    # CUDA.synchronize()
    return nothing
end
