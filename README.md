# Mapping Wrapper

A wrapper for batch processing the Mapping code, which takes command lines arugments to change either the density or the electron temperature.

## Example Command Line Arguments (in IPython)

Using default values T = 1e4, and density = 1:

```
run AutoMapping -numOfModels 50
```

Changing the density from 2 to 10 with constant temperature, T = 1e4:

```
run AutoMapping -numOfModels 50 -initDens 2  -finalDens 10
```

Changing the temperature from 1e5 to 1e6 with constant density = 1:

```
run AutoMapping -numOfModels 50 -initTemp 1e5  -finalTemp 1e6
```
