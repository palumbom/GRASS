# Installation

GRASS only requires a few steps to install. Simply clone the repo to your desired directory...

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
