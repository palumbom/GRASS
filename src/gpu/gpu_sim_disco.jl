function disk_sim_eclipse_disco_gpu(
    spec::SpecParams{T1}, disk::DiskParamsEclipse{T1}, 
    gpu_allocs::GPUAllocsEclipse{T2},
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

    μs = gpu_allocs.μs
    ld = gpu_allocs.ld
    ext = gpu_allocs.ext
    dA = gpu_allocs.dA
    z_rot = gpu_allocs.z_rot
    contrast = gpu_allocs.contrast

    n_patches = CUDA.length(μs)
    # DISCO supplies the line profile; LARS input profiles are not used on this path.
    # Preserve flux across templates; GPUAllocsEclipse initializes it to continuum.
    CUDA.fill!(prof, 0.0)

    # Thread/block configurations matching 2D CUDA grid conventions
    threads4 = (16, 16)
    blocks4 = (cld(n_patches, 16), cld(Nλ, 16))
    disco_centroid_offsets = CUDA.zeros(Float32, n_patches)

    threads5 = 1024
    blocks5 = cld(CUDA.length(prof), prod(threads5))
    threads_shift = 256
    blocks_shift = cld(n_patches, threads_shift)

    ext_toggle_val = ext_toggle_gpu ? T1(1) : T1(0)
    spot_toggle_val = spot_toggle_gpu ? T1(1) : T1(0)

    # Loop over time steps
    for t in 1:Nt
        # 1. Calculate eclipse quantities (patch μ, position, extinction mask, etc.)
        calc_eclipse_quantities_gpu!(
            time_stamps[t], obs_long, obs_lat, alt, wavelength, LD_type, 
            ext_toggle_val, ext_coeff, disk, gpu_allocs, spot_toggle_val
        )

        if skip_times[t]
            continue
        end

        # Seed per timestep for independent stochastic DISCO granulation draws per epoch
        epoch_seed = UInt32(0x12345678 + t)

        CUDA.@sync @cuda threads=threads_shift blocks=blocks_shift disco_centroid_offset_gpu!(
            disco_centroid_offsets, μs, disco_params, epoch_seed,
            Float32(spec.lines[1] / 10)
        )

        # Loop over spectral lines to synthesize
        for l in eachindex(spec.lines)
            if ext_toggle_gpu == true
                CUDA.@sync sum_wts = CUDA.sum(dA .* ld[:,:,l] .* ext[:,:,l])
            else
                CUDA.@sync sum_wts = CUDA.sum(dA .* ld[:,:,l])
            end
            # Reset line profile workspace to 0 before accumulating patch contributions
            CUDA.fill!(prof, 0.0)

            # 2. DISCO Line Profile Synthesis Kernel
            CUDA.@sync @cuda threads=threads4 blocks=blocks4 line_profile_disco_gpu!(
                l, prof, disco_centroid_offsets, μs, ld, dA, ext, λs, z_rot,
                contrast, ext_toggle_val, disco_params, epoch_seed
            )

            # 3. Normalize and accumulate spectrum frame into output flux matrix
            CUDA.@sync @cuda threads=threads5 blocks=blocks5 GRASS.apply_line!(t, prof, flux, sum_wts)
        end

    end

    # Copy output flux matrix from GPU to host CPU
    CUDA.@sync flux_cpu .= Array(flux)
    flux_cpu[:, skip_times] .= zero(eltype(flux_cpu))

    CUDA.synchronize()
    return nothing
end

function disco_centroid_offset_gpu!(
    offsets, μs, p_disco::DISCOParamsDevice, epoch_seed::UInt32, line_center
)
    idx = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    stride = blockDim().x * gridDim().x
    nθ_max = CUDA.size(μs, 2)
    n_wave = Int(p_disco.n_wave)
    wavelength = p_disco.wavelength

    for patch_id in idx:stride:CUDA.length(μs)
        row = (patch_id - 1) ÷ nθ_max
        col = (patch_id - 1) % nθ_max
        μ = μs[row + 1, col + 1]
        if μ <= 0.0f0
            @inbounds offsets[patch_id] = 0.0f0
            continue
        end

        μ_disco = Float32(μ)
        area = 0.0f0
        moment = 0.0f0
        f0 = disco_intensity(
            p_disco, μ_disco, Int32(1), Int32(patch_id), epoch_seed
        )
        for k in 1:(n_wave - 1)
            λ0 = wavelength[k]
            λ1 = wavelength[k + 1]
            f1 = disco_intensity(p_disco, μ_disco, Int32(k + 1), Int32(patch_id), epoch_seed)
            depth0 = max(1.0f0 - f0, 0.0f0)
            depth1 = max(1.0f0 - f1, 0.0f0)
            dλ = λ1 - λ0

            area += 0.5f0 * (depth0 + depth1) * dλ
            moment += 0.5f0 * (
                (λ0 - line_center) * depth0 + (λ1 - line_center) * depth1
            ) * dλ
            f0 = f1
        end
        @inbounds offsets[patch_id] = area > 0.0f0 ? moment / area : 0.0f0
    end

    return nothing
end

"""
    line_profile_disco_gpu!(...)

CUDA kernel for evaluating DISCO intensity per disk patch and accumulating into `prof`.
Applies local limb angle μ, rotation Doppler shift, limb darkening, and lunar occultation extinction.
Line shape and convective blueshift are naturally produced by DISCO.
"""
function line_profile_disco_gpu!(
    l, prof, centroid_offsets, μs, ld, dA, ext, λs, z_rot, contrast,
    ext_toggle, p_disco::DISCOParamsDevice, epoch_seed::UInt32
)
    idx = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    sdx = blockDim().x * gridDim().x
    idy = threadIdx().y + blockDim().y * (blockIdx().y - 1)
    sdy = blockDim().y * gridDim().y

    Nθ_max = CUDA.size(μs, 2)
    n_patches = CUDA.length(μs)
    n_λ = CUDA.length(λs)

    # Rest-frame limits and samples in nm from DISCO's wavelength table.
    rest_lo = p_disco.wavelength[1]
    rest_hi = p_disco.wavelength[Int(p_disco.n_wave)]

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
        μ_disco = Float32(μ)

        # Patch weight: dA * limb_darkening [* extinction_mask]
        w = dA[m, n] * ld[m, n, l]
        if ext_toggle == 1.0f0
            w *= ext[m, n, l]
        end

        w <= 0.0f0 && continue

        # Compute total Doppler shift factor (rotation)
        z_tot = (1.0f0 + z_rot[m, n])
        λ_spot_correction_nm =
            (1.0f0 - Float32(contrast[m, n])) * centroid_offsets[i]

        # Parallelized loop over observer wavelength grid points
        for j in idy:sdy:n_λ
            λ_obs = λs[j]
            # Scale the DISCO convective shift by the unspotted fraction.
            λ_rest_nm = Float32((λ_obs / 10.0) / z_tot) + λ_spot_correction_nm

            if λ_rest_nm >= rest_lo && λ_rest_nm <= rest_hi
                # DISCO's wavelength table is not uniform; locate the actual
                # bracketing samples instead of dividing by an average step.
                lo = Int32(1)
                hi = p_disco.n_wave
                while hi - lo > Int32(1)
                    mid = (lo + hi) ÷ Int32(2)
                    if p_disco.wavelength[Int(mid)] <= λ_rest_nm
                        lo = mid
                    else
                        hi = mid
                    end
                end
                t = (λ_rest_nm - p_disco.wavelength[Int(lo)]) /
                    (p_disco.wavelength[Int(lo + Int32(1))] -
                     p_disco.wavelength[Int(lo)])

                f0 = disco_intensity(p_disco, μ_disco, lo, Int32(i), epoch_seed)
                f1 = disco_intensity(
                    p_disco, μ_disco, lo + Int32(1), Int32(i), epoch_seed)

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