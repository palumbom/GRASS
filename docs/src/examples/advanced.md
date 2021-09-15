# Advanced Usage

Examples of more advanced usage of GRASS can be found in the [figures/ directory](https://github.com/palumbom/GRASS/tree/main/figures). These scripts were used to generate the figures that can be found in the GRASS paper. For clarity, a few specific advanced use cases are highlighted below.

## Accessing the Input Data
The "input data" described in the GRASS paper can be accessed as follows:

```julia
using GRASS

# get input data
input_data = GRASS.SolarData()
```

The  ```input_data``` variable is an instance of the ```SolarData``` composite type. It contains the bisector and width measurements for all epochs of the input data. The data are stored in dictionaries keyed by tuples of symbols representing a location on the disk. Accessing the data at a location on the disk is as simple as specifying the appropriate key.

```julia
# get bisector data at disk center
wav = input_data.wav[(:c, :mu10)]
bis = input_data.bis[(:c, :mu10)]
```

In the above code, ```wav``` is a 2D array of wavelengths. Each column contains a measurement of wavelength corresponding to the intensity measured in the corresponding column of ```bis``` at a given time. Plotting a bisector in matplotlib simply entails slicing out the desired column.

```julia
plt.plot(wav[:,1], bis[:,1])
plt.show()
```

## Measuring a Line Bisector
GRASS can also compute bisectors for absorption lines or CCF profiles. A detailed example of this can be seen in [the code for Figure 3](https://github.com/palumbom/GRASS/blob/main/figures/fig3.jl). A streamlined version is presented below. First, GRASS is used to generate some synthetic spectra (exactly as in the [basic use example](https://palumbom.github.io/GRASS/dev/examples/examples/#Generating-Synthetic-Spectra)):

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

Then, we call the ```GRASS.measure_bisector``` function:

```julia
wav_synth, bis_synth = GRASS.measure_bisector(wavelengths, flux, interpolate=true, top=0.9)
```

Two methods for calculating bisectors are implemented. By default, GRASS will attempt to interpolate the absorption line profile onto a uniform grid of flux values. Then the bisector is calculated as the average wavelength value on the left and right wings of the absorption line. Interpolation usually provides cleaner results, but it will fail if the line profile is sufficiently noisy. In this case, if ```interpolation=false```, GRASS will use an iterative approach that is more stable, but less precise.
