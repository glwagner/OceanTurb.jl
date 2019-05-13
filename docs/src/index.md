# OceanTurb.jl

> *How inappropriate to call this planet Earth when it is quite clearly Ocean.*
>> *Arthur C. Clark*

`OceanTurb.jl` implements one-dimensional partial differential equations
that model turbulent convection and diffusion in the ocean surface boundary layer.
It's purpose is the exploration, development, and practical usage
of ocean turbulence models.

### In the scheme of things

Just beneath the surface of the ocean, atmospheric fluxes of energy, heat
fresh water, salinity, and momentum
due to wind, waves, precipitation, evaporation, heating, cooling,
and radiation drive turbulence and mediate the exchange of quantities like
heat, momentum, and carbon between the atmosphere and ocean interior.
Models that approximate the effects of atmospheric forcing on
turbulence and turbulent mixing in the upper ocean are critical
components in ocean circulation models and coupled climate models.

## Installation

Open [julia](https://julialang.org), press `]` at
[the Julian prompt](https://docs.julialang.org/en/v1/stdlib/REPL/index.html#The-different-prompt-modes-1)
to enter
[package manager mode](https://docs.julialang.org/en/v1/stdlib/Pkg/#Pkg-1),
and type

```julia
pkg> add https://github.com/glwagner/OceanTurb.jl.git
```

This installs `OceanTurb.jl` from the master branch.
The tests can then be run by typing

```julia
pkg> test OceanTurb
```

Use help mode by typing `?` to find information about key functions:

```julia
help?> iterate!
```

which should give the information

```@docs
iterate!
```  

## Components

Solvers for various turbulence models are implemented in submodules
of `OceanTurb.jl`.
For example, our simplest model solves the 1D diffusion equation
and can be used by writing

```julia
julia> using OceanTurb.Diffusion
```

In addition to simple diffusion we have models for

* The K-Profile-Parameterization proposed by
* A 'modular' K-Profile-Parameterization that allows mixing and matching the three separate components of a canonical 'KPP'-type model.
* Pacanowski-Philander

## Authors

[Gregory L. Wagner](https://glwagner.github.io).
