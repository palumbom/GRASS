function disk_sim_eclipse_disco_gpu(
    spec::SpecParams{T1}, disk::DiskParamsEclipse{T1}, 
    soldata::GPUSolarData{T2}, gpu_allocs::GPUAllocsEclipse{T2},
    flux_cpu::AA{T1,2}, obs_long::T1, obs_lat::T1, alt::T1, time_stamps::Vector{String}, wavelength,
    ext_coeff, ext_toggle_gpu::Bool, spot_toggle_gpu::Bool, LD_type::String,
    disco_params::DISCOParamsDevice;
    skip_times::BitVector=falses(disk.Nt)
) where {T1<:AF, T2<:AF}

    # Dimensions for memory allocation
    Nt = disk.Nt
    Nλ = length(spec.lambdas)

    # Parse out GPU allocations
    λs = gpu_allocs.λs
    prof = gpu_allocs.prof
    flux = gpu_allocs.flux

    tloop = gpu_allocs.tloop
    dat_idx = gpu_allocs.dat_idx

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    contrast = gpu_allocs.contrast

    # Alias length array for random number loop bounds
    lenall_gpu = soldata.len

    n_patches = CUDA.length(μs)

    # Thread/block configurations matching 2D CUDA grid conventions
    threads1 = 1024
    blocks1 = cld(n_patches, prod(threads1))

    threads4 = (16, 16)
    blocks4 = (cld(n_patches, 16), cld(Nλ, 16))

    threads5 = 1024
    blocks5 = cld(CUDA.length(prof), prod(threads5))

    ext_toggle_val = ext_toggle_gpu ? 1.0 : 0.0
    spot_toggle_val = spot_toggle_gpu ? 1.0 : 0.0

    # Loop over time steps
    for t in 1:Nt
        # 1. Calculate eclipse quantities (patch μ, position, extinction mask, etc.)
        calc_eclipse_quantities_gpu!(
            time_stamps[t], obs_long, obs_lat, alt, wavelength, LD_type, 
            ext_toggle_val, ext_coeff, disk, gpu_allocs, spot_toggle_val
        )

        if isone(t)
            CUDA.@sync @cuda threads=threads1 blocks=blocks1 GRASS.generate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
        end

        if skip_times[t]
            CUDA.@sync @captured @cuda threads=threads1 blocks=blocks1 GRASS.iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
            continue
        end

        # Seed per timestep for independent stochastic DISCO granulation draws per epoch
        epoch_seed = UInt32(0x12345678 + t)

        # Loop over spectral lines to synthesize
        for l in eachindex(spec.lines)
            if ext_toggle_gpu == true
                CUDA.@sync sum_wts = CUDA.sum(dA .* ld[:,:,l] .* ext[:,:,l])
            else
                CUDA.@sync sum_wts = CUDA.sum(dA .* ld[:,:,l])
            end

            # 2. DISCO Line Profile Synthesis Kernel
            # Replaces trim_bisector_gpu!, fill_workspaces_2D_eclipse!, and line_profile_gpu!
            CUDA.@sync @cuda threads=threads4 blocks=blocks4 line_profile_disco_gpu!(
                l, prof, μs, ld, dA, ext, λs, z_rot, contrast,
                ext_toggle_val, disco_params, epoch_seed
            )

            # 3. Normalize and accumulate spectrum frame into output flux matrix
            CUDA.@sync @cuda threads=threads5 blocks=blocks5 GRASS.apply_line!(t, prof, flux, sum_wts)
        end

        # Iterate time loop index
        CUDA.@sync @captured @cuda threads=threads1 blocks=blocks1 GRASS.iterate_tloop_gpu!(tloop, dat_idx, lenall_gpu)
    end

    # Copy output flux matrix from GPU to host CPU
    CUDA.@sync flux_cpu .= Array(flux)

    CUDA.synchronize()
    return nothing
end

"""
    line_profile_disco_gpu!(...)

CUDA kernel for evaluating DISCO intensity per disk patch and accumulating into `prof`.
Applies local limb angle μ, rotation Doppler shift, limb darkening, and lunar occultation extinction.
Line shape and convective blueshift are naturally produced by DISCO.
"""
function line_profile_disco_gpu!(
    l, prof, μs, ld, dA, ext, λs, z_rot, contrast,
    ext_toggle, p_disco::DISCOParamsDevice, epoch_seed::UInt32
)
    idx = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y - 1)
    sdy = blockDim().y * gridDim().y

    Nθ_max = CUDA.size(μs, 2)
    n_patches = CUDA.length(μs)
    n_λ = CUDA.length(λs)

    # Rest frame limits from DISCO parameter model
    rest_lo = p_disco.wavelength[1]
    rest_hi = p_disco.wavelength[Int(p_disco.n_wave)]
    rest_step = (rest_hi - rest_lo) / Float32(p_disco.n_wave - Int32(1))

    # Parallelized loop over active disk patches
    for i in idx:sdx:n_patches
        row = (i - 1) ÷ Nθ_max
        col = (i - 1) % Nθ_max
        m = row + 1
        n = col + 1

        # Skip off-disk or occulted patches
        μ = μs[m, n]
        if μ <= 0.0f0
            continue
        end

        # Patch weight: dA * limb_darkening [* extinction_mask]
        w = dA[m, n] * ld[m, n, l]
        if ext_toggle == 1.0f0
            w *= ext[m, n, l]
        end

        w <= 0.0f0 && continue

        # Compute total Doppler shift factor (rotation)
        z_tot = (1.0f0 + z_rot[m, n])

        # Parallelized loop over observer wavelength grid points
        for j in idy:sdy:n_λ
            λ_obs = λs[j]
            # Convert observer wavelength to line rest frame
            λ_rest = λ_obs / z_tot

            pos = (λ_rest - rest_lo) / rest_step

            if pos >= 0.0f0 && pos <= Float32(p_disco.n_wave - Int32(1))
                k = Int32(floor(pos))
                k = min(k, p_disco.n_wave - Int32(2))
                t = pos - Float32(k)

                # Sample DISCO intensity at wavelength grid points k+1 and k+2
                f0 = disco_intensity(p_disco, Float32(μ), k + Int32(1), Int32(i), epoch_seed)
                f1 = disco_intensity(p_disco, Float32(μ), k + Int32(2), Int32(i), epoch_seed)

                # Linear interpolation in wavelength
                I_val = f0 + t * (f1 - f0)
            else
                # Outside line rest window: unabsorbed continuum (intensity = 1.0)
                I_val = 1.0f0
            end

            # Accumulate weighted contribution into global profile workspace
            CUDA.@atomic prof[j] += w * I_val
        end
    end

    return nothing
end