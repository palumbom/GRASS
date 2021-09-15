## Measuring Velocities from Synthetic Spectra

GRASS wraps the [EchelleCCFs package]() to measure apparent Doppler velocities from spectra.

```julia
# generate a velocity grid and compute the cross-correlation function
vel_grid, ccf = calc_ccf(wavelengths, flux, spec, normalize=true)

# measure radial velocities and uncertainties from the ccfs
rvs, sigs = calc_rvs_from_ccf(vel_grid, ccf)
```
