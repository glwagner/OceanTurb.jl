# OceanTurb.jl

`OceanTurb.jl` implements one-dimensional partial differential equations
that model turbulent convection and diffusion.

The primary purpose of the package is to explore models for convection
and turbulent mixing in the ocean's surface boundary layer, where atmospheric
forcing due to wind, waves, precipitation, evaporation, heating, cooling,
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
search: iterate! iterate InteractiveUtils isinteractive Iterators

  iterate!(model, Δt, nt=1)

  Step model forward in time by one time-step with step-size Δt.
```

## Components

Solvers for various turbulence models are implemented in submodules
of `OceanTurb.jl`.
For example, our simplest model solves the 1D diffusion equation
and can be used by writing

```julia
julia> using OceanTurb.Diffusion
```

In addition to simple diffusion we have

* K-Profile-Parameterization
* ... others coming soon.

## Authors

[Gregory L. Wagner](https://glwagner.github.io).
