# Public Functions 

```@meta
CurrentModule = GRASS
```

The public functions exported by GRASS are documented on this page. The high-level convenience functions  should meet the needs of most users. Some other potentially useful, lower-level methods are exposed by GRASS and documented below. A full index of all methods defined by GRASS is available in the [Full Index](@ref "Full Index").

## High-level Convenience Functions

GRASS provides a few high-level convenience wrappers for generating spectra. 

```@docs
synthesize_spectra
spec
disk
```


## Velocity Measurement 

GRASS wraps [EchelleCCFs.jl](https://github.com/RvSpectML/EchelleCCFs.jl) to measure velocities from spectra using the CCF method. 

```docs
calc_ccf
calc_rvs_from_ccf
```

## Input Data 

GRASS uses solar observations as input to simulate granulation in spectra. There are a few functions and composite types used to manipulate these data. See also [Input Data](@ref "Input Data").

```docs
LineProperties
SolarData
```

## Bisector Measurement

```docs
calc_bisector
calc_bisector_inverse_slope
calc_bisector_span
calc_bisector_bottom
calc_bisector_curvature
```

## Utilities 

```docs
calc_rms
```