# OceanTurb.jl

| **Documentation**             | **Build Status**                    | **License** |
|:-----------------------------:|:-----------------------------------:|:-----------:|
| [![docs][docs-img]][docs-url] | [![travis][travis-img]][travis-url] |[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://mit-license.org/)|


`OceanTurb.jl` provides software for solving one-dimensional 
models that approximate the physics of the 
ocean's turbulent surface boundary layer.

## Installation

Open [julia](https://julialang.org), press `]` to enter package manager mode, and type

```julia
pkg> add https://github.com/glwagner/OceanTurb.jl.git
```

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


[docs-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-url]: https://glwagner.github.io/OceanTurb.jl/latest/

[travis-img]: https://travis-ci.org/glwagner/OceanTurb.jl.svg?branch=master
[travis-url]: https://travis-ci.org/glwagner/OceanTurb.jl
