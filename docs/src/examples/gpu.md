# GPU Acceleration

GRASS ships with a GPU implementation of its disk synthesis that has been
validated to reproduce the fiducial CPU results within numerical precision. For
large simulations (many epochs, many lines, or fine disk grids) it can be
substantially faster than the CPU path.

!!! note
    The GPU implementation requires a functional CUDA device. [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl)
    is already a dependency of GRASS, so no additional installation is needed —
    but `synthesize_spectra` will throw an `AssertionError` (from
    `CUDA.functional()`) if no usable NVIDIA GPU is present.

## Enabling the GPU

The GPU is opt-in. Build `SpecParams` and `DiskParams` exactly as you would for a
CPU run (see [Basic Usage](@ref "Basic Usage")), then pass `use_gpu=true` to
[`synthesize_spectra`](@ref):

```julia
using GRASS

# the same parameters used for a CPU synthesis
λrest = 5434.5
spec = SpecParams(lines=[λrest], depths=[0.75], resolution=7e5)
disk = DiskParams(Nt=50)

# run the synthesis on the GPU
wavelengths, flux = synthesize_spectra(spec, disk, use_gpu=true)
```

The return values and their shapes are identical to the CPU path, so every
downstream measurement — velocities, bisectors, and other diagnostics — works
without modification.

## Numerical Precision

By default the GPU implementation uses **double-precision** (`Float64`) floats.
This is the safe and recommended choice. Single precision is supported via the
`precision` keyword, but it is **not recommended**:

```julia
# NOT recommended — single precision triggers a warning and large errors
wavelengths, flux = synthesize_spectra(spec, disk, use_gpu=true, precision=Float32)
```

Requesting `Float32` emits a warning at runtime:

```
┌ Warning: Single-precision GPU implementation produces large flux and velocity errors!
```

This is because catastrophic cancellation in an internal interpolation operation
produces large flux errors in single precision. See the [Caveats](@ref "Caveats")
page for the quantitative impact. Double precision avoids the problem entirely,
at the cost of some performance on hardware with limited double-precision
throughput (e.g. many consumer GPUs).

## Eclipse Synthesis on the GPU

The eclipse entry point [`synthesize_spectra_eclipse`](@ref) accepts the same
`use_gpu` and `precision` keywords, with the same precision guidance:

```julia
# spec, disk (a DiskParamsEclipse), and the observer arguments are built exactly
# as in the Eclipse Synthesis tutorial; only use_gpu is added here.
wavelengths, flux = synthesize_spectra_eclipse(spec, disk, [λrest], LD_type,
                                               obs_long, obs_lat, alt,
                                               time_stamps, ext_coeff,
                                               use_gpu=true)
```

See [Eclipse Synthesis](@ref "Eclipse Synthesis") for the full eclipse workflow.
