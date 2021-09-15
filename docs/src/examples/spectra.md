The simplest use case for GRASS is the generation of time series of synthetic spectra. This can be done in only a few lines of code. The following example generates 25 spectra consisting of a single line at 5434.5 Angstroms. By default, the temporal spacing of the spectra is 15 seconds.

```julia
using GRASS

# parameters for lines in the spectra
lines = [5434.5]
depths = [0.75]
resolution = 7e5
spec = SpecParams(lines=lines, depths=depths, resolution=resolution)

# specify number of epochs (default 15-second spacing)
disk = DiskParams(Nt=25)

# synthesize the spectra
wavelengths, flux = synthesize_spectra(spec, disk)
```
