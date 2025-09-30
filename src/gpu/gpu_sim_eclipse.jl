function disk_sim_eclipse_gpu(spec::SpecParams{T1}, disk::DiskParamsEclipse{T1}, 
                               soldata::GPUSolarData{T2}, gpu_allocs::GPUAllocsEclipse{T2},
                               flux_cpu::AA{T1,2}, 
                               obs_long::T1, obs_lat::T1, alt::T1, time_stamps::Vector{Float64}, wavelength,
                               ext_coeff, ext_toggle_cpu::Bool, LD_type::String;
                               skip_times::BitVector=falses(disk.Nt)) where {T1<:AF, T2<:AF}

    # get dimensions for memory alloc
    Nt = disk.Nt
    Nλ = length(spec.lambdas)                           
                                
    # parse out GPU allocations
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux

    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints
    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs

    # alias the input data from GPUSolarData
    lenall_gpu = soldata.len
    bisall_gpu = soldata.bis
    intall_gpu = soldata.int
    widall_gpu = soldata.wid
    depcontrast_gpu = soldata.dep_contrast

    # set number of threads and blocks for len(μ) gpu kernels
    threads1 = 1024
    blocks1 = cld(CUDA.length(μs), prod(threads1))

    # set number of threads and blocks for trimming functions
    threads2 = (4,4,16)
    blocks2 = cld(length(lenall_gpu) * maximum(lenall_gpu) * 100, prod(threads2))

    # set number of threads and blocks for len(μ) * 100 matrix gpu functions
    threads3 = (16,16)
    blocks3 = cld(CUDA.length(μs) * 100, prod(threads3))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads4 = (16,16)
    blocks4 = cld(CUDA.length(μs) * Nλ, prod(threads4))

    threads5 = 1024
    blocks5 = cld(CUDA.length(prof), prod(threads5))

    # allocate destinations for interpolations
    @cusync begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
    end

    # allocate memory for means
    @cusync begin
        bisall_mean = CUDA.zeros(CUDA.eltype(bisall_gpu_loop), 100, CUDA.size(bisall_gpu_loop, 3))
        intall_mean = CUDA.zeros(CUDA.eltype(intall_gpu_loop), 100, CUDA.size(intall_gpu_loop, 3))
        widall_mean = CUDA.zeros(CUDA.eltype(widall_gpu_loop), 100, CUDA.size(widall_gpu_loop, 3))
    end

    threads6 = (4, 16)
    blocks6 = cld(length(lenall_gpu) * 100, prod(threads6))

    @cusync @cuda threads=threads6 blocks=blocks6 time_average_bis!(lenall_gpu, bisall_mean, intall_mean, 
                                                                    widall_mean, bisall_gpu, intall_gpu, 
                                                                    widall_gpu)

    if ext_toggle_cpu == true
        ext_toggle_gpu = 1.0
    else 
        ext_toggle_gpu = 0.0
    end

    brightness = Vector{Float64}(undef,size(time_stamps)...)
    # loop over time
    for t in 1:Nt
        # sort out the system geometry
        if ext_coeff == "on"
            if t < 48
                ext_coeff1 = ext_coeff_file[ext_coeff_file[!, "Wavelength"] .== wavelength, "ext1"]
                # compute intensity for timestamp
                calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, wavelength, LD_type, ext_toggle_gpu, ext_coeff1[1], disk, gpu_allocs)
            elseif t >= 48
                ext_coeff2 = ext_coeff_file[ext_coeff_file[!, "Wavelength"] .== wavelength, "ext2"]
                # compute intensity for timestamp
                calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, wavelength, LD_type, ext_toggle_gpu, ext_coeff2[1], disk, gpu_allocs)
            end
        else
            calc_eclipse_quantities_gpu!(time_stamps[t], obs_long, obs_lat, alt, wavelength, LD_type, ext_toggle_gpu, ext_coeff, disk, gpu_allocs)
        end

        # get conv. blueshift and keys from input data
        get_keys_and_cbs_gpu!(gpu_allocs, soldata)

        if isone(t)
            # generate the random numbers on the gpu
            @cusync @cuda threads=threads1 blocks=blocks1 generate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
        end

        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            if ext_toggle_cpu == true
                # get weighted disk average cbs
                @cusync sum_wts = CUDA.sum(dA .* ld[:,:,l] .* ext[:,:,l])
                @cusync z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l] .* ext[:,:,l]) / sum_wts
                brightness[t] = sum_wts
            else
                # get weighted disk average cbs
                @cusync sum_wts = CUDA.sum(dA .* ld[:,:,l])
                @cusync z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l]) / sum_wts
                brightness[t] = sum_wts
            end

            # calculate how much extra shift is needed
            extra_z = spec.conv_blueshifts .- z_cbs_avg

            # trim all the bisector data
            @cusync @cuda threads=threads2 blocks=blocks2 trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                             depcontrast_gpu, lenall_gpu,
                                                                             bisall_gpu_loop, intall_gpu_loop,
                                                                             widall_gpu_loop, bisall_gpu,
                                                                             intall_gpu, widall_gpu)

            # assemble line shape on even int grid
            @cusync @cuda threads=threads3 blocks=blocks3 fill_workspaces_2D_eclipse!(spec.lines[l], spec.variability[l], extra_z[l],
                                                                           tloop, dat_idx,
                                                                           z_rot, z_cbs, lenall_gpu,
                                                                           bisall_gpu_loop, intall_gpu_loop,
                                                                           widall_gpu_loop, allwavs, allints)
            
            # do the line synthesis, interp back onto wavelength grid
            @cusync @cuda threads=threads4 blocks=blocks4 line_profile_gpu!(l, prof, μs, ld, dA, ext, λs, allwavs, allints, ext_toggle_gpu)

            # copy data from GPU to CPU
            @cusync @cuda threads=threads5 blocks=blocks5 apply_line!(t, prof, flux, sum_wts)
        end

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # copy over flux
    @cusync flux_cpu .= Array(flux)

    # make sure nothing is still running on GPU
    CUDA.synchronize()
    return
end

function disk_sim_eclipse_gpu(spec::SpecParams{T1}, disk::DiskParamsEclipse{T1}, 
                               soldata::GPUSolarData{T2}, gpu_allocs::GPUAllocsEclipse{T2},
                               flux_cpu::AA{T1,2}, 
                               obs_long::T1, obs_lat::T1, alt::T1, time_stamps::Vector{Float64}, wavelength,
                               ext_coeff, ext_toggle_cpu::Bool, LD_type::String, LD1, LD2, LD3, LD4, CB1, CB2, CB3;
                               skip_times::BitVector=falses(disk.Nt)) where {T1<:AF, T2<:AF}

    # get dimensions for memory alloc
    Nt = disk.Nt
    Nλ = length(spec.lambdas)                           
                                
    # parse out GPU allocations
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux

    allwavs = gpu_allocs.allwavs
    allints = gpu_allocs.allints
    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    z_cbs = gpu_allocs.z_cbs

    # alias the input data from GPUSolarData
    lenall_gpu = soldata.len
    bisall_gpu = soldata.bis
    intall_gpu = soldata.int
    widall_gpu = soldata.wid
    depcontrast_gpu = soldata.dep_contrast

    # set number of threads and blocks for len(μ) gpu kernels
    threads1 = 1024
    blocks1 = cld(CUDA.length(μs), prod(threads1))

    # set number of threads and blocks for trimming functions
    threads2 = (4,4,16)
    blocks2 = cld(length(lenall_gpu) * maximum(lenall_gpu) * 100, prod(threads2))

    # set number of threads and blocks for len(μ) * 100 matrix gpu functions
    threads3 = (16,16)
    blocks3 = cld(CUDA.length(μs) * 100, prod(threads3))

    # set number of threads and blocks for N*N*Nλ matrix gpu functions
    threads4 = (16,16)
    blocks4 = cld(CUDA.length(μs) * Nλ, prod(threads4))

    threads5 = 1024
    blocks5 = cld(CUDA.length(prof), prod(threads5))

    # allocate destinations for interpolations
    @cusync begin
        bisall_gpu_loop = CUDA.copy(bisall_gpu)
        intall_gpu_loop = CUDA.copy(intall_gpu)
        widall_gpu_loop = CUDA.copy(widall_gpu)
    end

    # allocate memory for means
    @cusync begin
        bisall_mean = CUDA.zeros(CUDA.eltype(bisall_gpu_loop), 100, CUDA.size(bisall_gpu_loop, 3))
        intall_mean = CUDA.zeros(CUDA.eltype(intall_gpu_loop), 100, CUDA.size(intall_gpu_loop, 3))
        widall_mean = CUDA.zeros(CUDA.eltype(widall_gpu_loop), 100, CUDA.size(widall_gpu_loop, 3))
    end

    threads6 = (4, 16)
    blocks6 = cld(length(lenall_gpu) * 100, prod(threads6))

    @cusync @cuda threads=threads6 blocks=blocks6 time_average_bis!(lenall_gpu, bisall_mean, intall_mean, 
                                                                    widall_mean, bisall_gpu, intall_gpu, 
                                                                    widall_gpu)

    if ext_toggle_cpu == true
        ext_toggle_gpu = 1.0
    else 
        ext_toggle_gpu = 0.0
    end

    brightness = Vector{Float64}(undef,size(time_stamps)...)
    # loop over time
    for t in 1:Nt
        # sort out the system geometry
        if ext_coeff == "on"
            if t < 48
                ext_coeff1 = ext_coeff_file[ext_coeff_file[!, "Wavelength"] .== wavelength, "ext1"]
                # compute intensity for timestamp
                calc_eclipse_quantities_gpu!(time_stamps[t], t, obs_long, obs_lat, alt, wavelength, ext_toggle_gpu, ext_coeff1[1], disk, gpu_allocs, LD1, LD2, LD3, LD4, CB1, CB2, CB3)
            elseif t >= 48
                ext_coeff2 = ext_coeff_file[ext_coeff_file[!, "Wavelength"] .== wavelength, "ext2"]
                # compute intensity for timestamp
                calc_eclipse_quantities_gpu!(time_stamps[t], t, obs_long, obs_lat, alt, wavelength, ext_toggle_gpu, ext_coeff2[1], disk, gpu_allocs, LD1, LD2, LD3, LD4, CB1, CB2, CB3)
            end
        else
            calc_eclipse_quantities_gpu!(time_stamps[t], t, obs_long, obs_lat, alt, wavelength, ext_toggle_gpu, ext_coeff, disk, gpu_allocs, LD1, LD2, LD3, LD4, CB1, CB2, CB3)
        end

        # get conv. blueshift and keys from input data
        get_keys_and_cbs_gpu!(gpu_allocs, soldata)

        if isone(t)
            # generate the random numbers on the gpu
            @cusync @cuda threads=threads1 blocks=blocks1 generate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
        end

        # don't synthesize spectrum if skip_times is true, but iterate t index
        if skip_times[t]
            @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
        end

        # loop over lines to synthesize
        for l in eachindex(spec.lines)
            if ext_toggle_cpu == true
                # get weighted disk average cbs
                @cusync sum_wts = CUDA.sum(dA .* ld[:,:,l] .* ext[:,:,l])
                @cusync z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l] .* ext[:,:,l]) / sum_wts
                brightness[t] = sum_wts
            else
                # get weighted disk average cbs
                @cusync sum_wts = CUDA.sum(dA .* ld[:,:,l])
                @cusync z_cbs_avg = CUDA.sum(z_cbs .* dA .* ld[:,:,l]) / sum_wts
                brightness[t] = sum_wts
            end

            # calculate how much extra shift is needed
            extra_z = spec.conv_blueshifts .- z_cbs_avg

            # trim all the bisector data
            @cusync @cuda threads=threads2 blocks=blocks2 trim_bisector_gpu!(spec.depths[l], spec.variability[l],
                                                                             depcontrast_gpu, lenall_gpu,
                                                                             bisall_gpu_loop, intall_gpu_loop,
                                                                             widall_gpu_loop, bisall_gpu,
                                                                             intall_gpu, widall_gpu)

            # assemble line shape on even int grid
            @cusync @cuda threads=threads3 blocks=blocks3 fill_workspaces_2D_eclipse!(spec.lines[l], spec.variability[l], extra_z[l],
                                                                           tloop, dat_idx,
                                                                           z_rot, z_cbs, lenall_gpu,
                                                                           bisall_gpu_loop, intall_gpu_loop,
                                                                           widall_gpu_loop, allwavs, allints)
            
            # do the line synthesis, interp back onto wavelength grid
            @cusync @cuda threads=threads4 blocks=blocks4 line_profile_gpu!(l, prof, μs, ld, dA, ext, λs, allwavs, allints, ext_toggle_gpu)

            # copy data from GPU to CPU
            @cusync @cuda threads=threads5 blocks=blocks5 apply_line!(t, prof, flux, sum_wts)
        end

        # iterate tloop
        @cusync @captured @cuda threads=threads1 blocks=blocks1 iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # copy over flux
    @cusync flux_cpu .= Array(flux)

    # make sure nothing is still running on GPU
    CUDA.synchronize()
    return
end
