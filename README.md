# GRASS - GRanulation and Spectrum Simulator 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://palumbom.github.io/GRASS/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://palumbom.github.io/GRASS/dev)
[![Build Status](https://github.com/palumbom/GRASS/workflows/CI/badge.svg)](https://github.com/palumbom/GRASS/actions)
[![arXiv](https://img.shields.io/badge/arXiv-2110.11839-b31b1b.svg)](https://arxiv.org/abs/2110.11839)
[![arXiv](https://img.shields.io/badge/arXiv-2405.07945-b31b1b.svg)](https://arxiv.org/abs/2405.07945)

GRASS is a Juliavpackage designed to produce time series of stellar spectra with realistic line-shape changes from solar granulation. GRASS v1.0.x is described in detail in [Palumbo et al. (2022)](https://arxiv.org/abs/2110.11839); the results of this paper can be reproduced using the [showyourwork workflow](https://github.com/showyourwork/showyourwork) from [this repo](https://github.com/palumbom/palumbo22). GRASS v2.0.x is presented in [Palumbo et al. (2024a)](https://arxiv.org/abs/2405.07945); [this repo](https://github.com/palumbom/grass-two) reproduces the figures and quantitative results presented therein. 

## Installation

GRASS is written entirely in Julia and requires Julia v1.9 or greater. Installation instructions for Julia are available from [julialang.org](https://julialang.org/downloads/).

GRASS itself only requires a few steps to install. Simply clone the repo to your desired directory...

```bash
git clone git@github.com:palumbom/GRASS.git
cd GRASS 
julia
```

... and then add it with Julia's built-in package manager, `Pkg`:

```julia
using Pkg
Pkg.add(path=".") # assuming you are in /PATH/TO/GRASS
using GRASS
```

If you wish to develop or otherwise contribute to GRASS, instead add the package in develop mode:

```julia
using Pkg
Pkg.develop(path=".") # assuming you are in /PATH/TO/GRASS
using GRASS
```

Upon first invocation of GRASS, Julia will automatically install the package dependencies and download the required input data. The input data can be re-installed by invoking

```julia
Pkg.build("GRASS")
```

Alternatively, these data can be directly downloaded from [Zenodo](https://zenodo.org/records/8271417).

## Basic Example
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
[![DOI](https://zenodo.org/badge/364662564.svg)](https://zenodo.org/badge/latestdoi/364662564)
<a href="https://ascl.net/2110.011"><img src="https://img.shields.io/badge/ascl-2110.011-blue.svg?colorB=262255" alt="ascl:2110.011" /></a>
[![arXiv](https://img.shields.io/badge/arXiv-2110.11839-b31b1b.svg)](https://arxiv.org/abs/2110.11839)
[![arXiv](https://img.shields.io/badge/arXiv-2405.07945-b31b1b.svg)](https://arxiv.org/abs/2405.07945)

If you use GRASS in your research, please cite the relevant [software release](https://zenodo.org/badge/latestdoi/364662564) and [paper(s)](https://arxiv.org/abs/2110.11839).

## Author & Contact 
[![GitHub followers](https://img.shields.io/github/followers/palumbom?label=Follow&style=social)](https://github.com/palumbom)

This repo is maintained by Michael Palumbo. You may may contact him via his email - [palumbo@psu.edu](mailto:palumbo@psu.edu)
