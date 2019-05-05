# Mapping Wrapper

A wrapper for batch processing the Mapping code, which takes command lines arugments to change either the density or the electron temperature.

## Example Command Line Arguments (in IPython)

Using default values T = 1e4, and density = 1, for 50 models (which doesn't make sense because they would be the same):

```
run AutoMapping -numOfModels 50
```

Changing the density from 2 to 10 with constant temperature, T = 1e4:

```
run AutoMapping -numOfModels 500 -initDens 2  -finalDens 10
```

The increment to update the density is calculated, (finalDens - initDens) / (numOfModels-1), hence the increment for the density in the above example is 0.016.


Changing the temperature from 1e5 to 1e6 with constant density = 1:

```
run AutoMapping -numOfModels 500 -initTemp 1e5  -finalTemp 1e6
```

The increment to update the density is calculated, (finalTemp - initTemp) / (numOfModels-1), hence the increment for the density in the above example is ~ 2000[K].

Changing the temperature from 1e5 to 1e6 with density changing from 2 to 10:

```
run AutoMapping -numOfModels 500 -initTemp 1e5  -finalTemp 1e6  -initDens 2  -finalDens 10
```

The increment to update the density is calculated, (finalTemp - initTemp) / (numOfModels-1), hence the increment for the density in the above example is ~ 2000[K].


This code currently only changes the temperature and the density from the mapping defaults. All sources are switched off, (e.g. cosmic radiation heating), or set to the default behaviours for the code.
