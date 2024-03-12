# Basic Usage

The most basic use case of GRASS is generating time series of synthetic spectra and measuring apparent velocities. These can both be done in only a few lines of code.

## Generating Synthetic Spectra

The simplest use case for GRASS is the generation of time series of synthetic spectra. This can be done in only a few lines of code. The following example generates 25 spectra consisting of a single line at 5434.5 Angstroms. By default, the temporal spacing of the spectra is 15 seconds.

```julia
using GRASS

# parameters for lines in the spectra
lines = [5434.5]     # array of line centers in angstroms
depths = [0.75]      # continuum-normalized depth of lines
resolution = 7e5     # spectral resolution of the output spectra
spec = SpecParams(lines=lines, depths=depths, resolution=resolution)

# specify number of epochs (default 15-second spacing)
disk = DiskParams(Nt=25)

# synthesize the spectra
wavelengths, flux = synthesize_spectra(spec, disk)
```

## Measuring Velocities

GRASS wraps the [EchelleCCFs package](https://github.com/RvSpectML/EchelleCCFs.jl) to measure apparent Doppler velocities from spectra.

```julia
# generate a velocity grid and compute the cross-correlation function
vel_grid, ccf = calc_ccf(wavelengths, flux, spec, normalize=true)

# measure radial velocities and uncertainties from the ccfs
rvs, sigs = calc_rvs_from_ccf(vel_grid, ccf)
```
