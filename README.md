# OceanTurb.jl

This code provides models that approximate the physics of the 
turbulent surface boundary layer of the ocean.

## Installation

Open julia, press `]` to enter package manager mode, and type

```julia
pkg> add https://github.com/glwagner/OceanTurb.jl.git
```

## Documentation

* [**latest**](https://glwagner.github.io/OceanTurb.jl/latest)

## Example(s)

Check our 
[diffusion examples notebook](https://github.com/glwagner/OceanTurb.jl/blob/master/examples/diffusion_example.ipynb), 
or for the more hardcore, 
[our notebook of KPP examples](https://github.com/glwagner/OceanTurb.jl/blob/master/examples/kpp_examples.ipynb).
Even more example scripts and notebooks are eagerly waiting in `examples/`.

# The turbulence models

This code aims to provide solvers for as many ocean surface boundary layer
turbulence models as possible.
Currently we have solvers for

* the diffusion equation
* the [K-Profile-Parameterization](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
* [Pacanowski-Philander](https://journals.ametsoc.org/doi/abs/10.1175/1520-0485(1981)011%3C1443:POVMIN%3E2.0.CO;2)

# Authors

[Gregory Wagner](glwagner.github.io).
