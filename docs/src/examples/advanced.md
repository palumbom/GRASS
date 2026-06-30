# Manipulating Input Data

Examples of more advanced usage of GRASS can be found in the companion repositories that reproduce the figures in the GRASS papers: [palumbom/palumbo22](https://github.com/palumbom/palumbo22) for [Palumbo et al. (2022)](https://arxiv.org/abs/2110.11839) and [palumbom/grass-two](https://github.com/palumbom/grass-two) for [Palumbo et al. (2024a)](https://arxiv.org/abs/2405.07945). For clarity, a few specific advanced use cases are highlighted below.

## Synthesizing Multiple Lines

[`SpecParams`](@ref) accepts arrays, so synthesizing several lines at once is just a matter of passing one entry per line. Each line needs a center (in air wavelength, angstroms) and a continuum-normalized depth between 0 and 1:

```julia
using GRASS

lines  = [5432.5, 5434.5]   # line centers (Å)
depths = [0.40,   0.75]     # continuum-normalized depths
spec   = SpecParams(lines=lines, depths=depths, resolution=7e5)
disk   = DiskParams(Nt=50)

wavelengths, flux = synthesize_spectra(spec, disk)
```

The output wavelength grid automatically spans all requested lines (plus a small `buffer`, controllable via the `buffer` keyword to `SpecParams`).

## Choosing Input Templates

GRASS models each synthetic line by replaying the observed variability of a real solar line ("template"). If you do not specify a template, GRASS picks one automatically by finding the observed line closest to your requested line in depth and Landé g-factor — and it emits a warning to let you know it made that choice for you:

```
┌ Warning: No line template specified!
```

To choose templates explicitly, pass the `templates` keyword with one `.h5` filename per line. The available template files can be listed via [`LineProperties`](@ref):

```julia
# list the template lines bundled with GRASS
lp = GRASS.LineProperties()
GRASS.get_file(lp)     # template .h5 filenames
GRASS.get_depth(lp)    # their intrinsic depths
GRASS.get_geff(lp)     # their Landé g-factors
```

```julia
# synthesize the 5434.5 Å line using a specific template
spec = SpecParams(lines=[5434.5], depths=[0.75], templates=["FeI_5434"])
```

!!! warning
    A template observed for one line is not guaranteed to faithfully reproduce
    variability in an arbitrary other line, especially far from its observed
    wavelength or depth. See the [Caveats](@ref "Caveats") page before relying on
    cross-line templates.

## Toggling Variability per Line

By default every line carries granulation-driven variability. The per-line `variability` flag lets you freeze individual lines — useful for isolating the contribution of a single line, or for generating a static reference spectrum. Pass one Boolean per line:

```julia
# make the first line static and the second variable
spec = SpecParams(lines=[5432.5, 5434.5], depths=[0.40, 0.75],
                  variability=[false, true])
```

A static line still has the correct mean shape and convective blueshift; it simply does not vary in time.

## Accessing the Input Data

The "input data" described in the GRASS papers can be accessed directly through the [`SolarData`](@ref) composite type:

```julia
using GRASS

# get input data
input_data = GRASS.SolarData()
```

The `input_data` variable contains the bisector and width measurements for all epochs of the input data. The data are stored in dictionaries keyed by a tuple of two `Symbol`s: the disk axis (e.g. `:c` for center, plus the cardinal directions) and the center-to-limb position `mu` (e.g. `:mu10` for μ = 1.0 at disk center, down to the limb). The relevant fields are `bis` (bisector wavelengths) and `int` (the corresponding intensities); `wid` holds the line widths. Accessing the data at a location on the disk is as simple as specifying the appropriate key:

```julia
# get bisector data at disk center
bis = input_data.bis[(:c, :mu10)]   # 2D array of bisector wavelengths
int = input_data.int[(:c, :mu10)]   # 2D array of corresponding intensities
```

Both `bis` and `int` are 2D arrays whose columns are individual epochs: column `i` of `bis` gives the bisector wavelengths at the intensity levels in column `i` of `int`. Plotting a single bisector with matplotlib simply entails slicing out the desired column:

```julia
plt.plot(bis[:,1], int[:,1])
plt.show()
```
