# Manipulating Input Data

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