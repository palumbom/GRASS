# GRASS

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://palumbom.github.io/GRASS/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://palumbom.github.io/GRASS/dev)
[![Build Status](https://github.com/palumbom/GRASS/workflows/CI/badge.svg)](https://github.com/palumbom/GRASS/actions)
[![DOI](https://zenodo.org/badge/364662564.svg)](https://zenodo.org/badge/latestdoi/364662564)

GRASS is a package designed to produce realistic time series of stellar spectra with realistic line-shape changes from solar-like granulation.

## Introduction

## Examples
Generating synthetic spectra with GRASS only takes a few lines of Julia:

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

Additional details and examples can be found in [the documentation](https://palumbom.github.io/GRASS/stable).

## Citation

If you use GRASS in your research, please cite the relevant software release and paper(s).

## Author & Contact [![Twitter Follow](https://img.shields.io/twitter/follow/michael_palumbo?style=social)](https://twitter.com/michael_palumbo) [![GitHub followers](https://img.shields.io/github/followers/palumbom?label=Follow&style=social)](https://github.com/palumbom)

This repo is maintained by Michael Palumbo. You may may contact him via his email - palumbo@psu.edu
