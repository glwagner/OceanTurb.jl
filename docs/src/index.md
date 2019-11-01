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
pkg> add OceanTurb
```

Use help mode by typing `?` to find information about key functions:

```julia
help?> iterate!
```

which should give the information

```@docs
iterate!
```  

## Modules and models

Solvers for various turbulence models are implemented in submodules
of `OceanTurb.jl`.
For example, our simplest module solves the 1D diffusion equation.
A diffusion `Model` is instantiated by writing

```julia
using OceanTurb

# 100-grid point model with height 1.0 and diffusivity 0.01.
model = Diffusion.Model(grid = UniformGrid(N=100, L=1.0), parameters = Parameters(K=0.01))
```

Setting an initial condition is done by writing

```julia
c₀(z) = exp(-(z + 0.5)^2 / 0.005)
model.solution.c = c₀
```

Time stepping a model forward looks like

```julia
iterate!(model, Δt=0.01, Nt=100)
```

This example, and more, can be found in the `/examples` directory.

In addition to simple diffusion we have models for

* The K-Profile-Parameterization proposed by
    [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)

* A 'modular', and therefore generic, ``K``-profile parameterization that has multiple models
    for diffusivity, diffusivity shapes and profiles, nonlocal fluxes including a diagnostic plume model, 
    and mixing depth.

* The Pacanowski-Philander parameterization

## Authors

The author of this software is [Gregory L. Wagner](https://glwagner.github.io).
