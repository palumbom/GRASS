# Installation

GRASS is written entirely in Julia and requires Julia v1.9 or greater. Installation instructions for Julia are available from [julialang.org](https://julialang.org/downloads/); alternatively, Julia can be trivially installed on a Mac via [homebrew](https://formulae.brew.sh/cask/julia).

The installation of GRASS itself only requires a few steps to install. Simply clone the repo to your desired directory...

```bash
cd DIRECTORY
git clone git@github.com:palumbom/GRASS.git
```

... and then add it with Julia's Pkg:

```julia
using Pkg
Pkg.add(path="DIRECTORY")
using GRASS
```

If you wish to develop or otherwise contribute to GRASS, instead add the package in develop mode:

```julia
using Pkg
Pkg.develop(path="DIRECTORY")
using GRASS
```

Upon first invocation of GRASS, Julia will automatically install the package dependencies and download the required input data. The input data can be re-installed by invoking

```julia
Pkg.build("GRASS")
```

Alternatively, these data can be directly downloaded from [Zenodo](https://zenodo.org/records/8271417).
