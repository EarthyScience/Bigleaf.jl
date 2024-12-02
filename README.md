# Bigleaf

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://earthyscience.github.io/Bigleaf.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://earthyscience.github.io/Bigleaf.jl/dev/)
[![Build Status](https://github.com/EarthyScience/Bigleaf.jl/workflows/CI/badge.svg)](https://github.com/EarthyScience/Bigleaf.jl/actions)
[![Coverage](https://codecov.io/gh/bgctw/Bigleaf.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bgctw/Bigleaf.jl)


**Bigleaf.jl** is a partial Julia port of Jürgen Knauer's 
[**bigleaf** R package](https://bitbucket.org/juergenknauer/Bigleaf) 
for the calculation of physical (e.g. aerodynamic conductance, surface temperature) 
and physiological (e.g. canopy conductance, water-use efficiency) ecosystem properties 
from eddy covariance data and accompanying meteorological measurements. 
All calculations are based on a 'big-leaf' representation of the vegetation and 
return representative bulk ecosystem/canopy variables.


## Citation
Knauer J, El-Madany TS, Zaehle S, Migliavacca M (2018) Bigleaf—An R package for the calculation of physical and physiological ecosystem properties from eddy covariance data.
PLoS ONE 13(8): e0201114. [https://doi.org/10.1371/journal.pone.0201114](https://doi.org/10.1371/journal.pone.0201114)


## Installation and Loading

The `Bigleaf.jl` package can be installed with the usual command once:

```julia
using Pkg
Pkg.add("Bigleaf")
```

And then imported to every Julia session by:
```julia
using Bigleaf
```

## Usage
See the [Documentation](https://earthyscience.github.io/Bigleaf.jl/dev/) 
that includes a [walkthrough](https://earthyscience.github.io/Bigleaf.jl/dev/walkthrough/).

[Please report bugs or issues here](https://github.com/EarthyScience/Bigleaf.jl/issues)

## Package content 
We are porting functionality of the [R package](https://bitbucket.org/juergenknauer/Bigleaf) as needed. Please
file an [issue](https://github.com/EarthyScience/Bigleaf.jl/issues) if you need a specific feature.

At the current state we ported
- Meteorological variables
- Boundary layer and Aerodynamic conductance
- Surface conductance
- Evapotranspiration
- Potential radiation
- Unit conversions
- Filtering data

