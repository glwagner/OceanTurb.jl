```math
\newcommand{\c}{\, ,}
\newcommand{\p}{\, .}
\newcommand{\d}{\partial}

\newcommand{\r}[1]{\mathrm{#1}}

\newcommand{\ee}{\mathrm{e}}

\newcommand{\beq}{\begin{equation}}
\newcommand{\eeq}{\end{equation}}

\newcommand{\beqs}{\begin{gather}}
\newcommand{\eeqs}{\end{gather}}
```

# Physics and modeling of the oceanic boundary layer

Models for the oceanic boundary layer are partial
differential equations that approximate the effects of 

1. internal and surface fluxes, and
2. turbulent vertical fluxes

on the evolution of temperature, salinity, and momentum in
the near-surface ocean.

Internal and surface fluxes of heat are due to

* absorption of incoming shortwave solar radiation;
* cooling by outgoing longwave radiation;
* latent and sensible heat exchange with the atmosphere.

Surface fluxes of salinity occur due to evaporation and precipitation, 
while momentum fluxes are associated with atmospheric winds.

Vertical turbulent fluxes are typically associated with

* gravitational instability and convection, and
* mechanical turbulent mixing associated with currents and wind forcing.

`OceanMixedLayerModels.jl` uses 
an implementation of atmospheric and radiative forcings that is shared
across all models. The models therefore differ mainly in the way they
parameterize convective and mechanical mixing.

## Coordinate system

We use a Cartesian coordinate system in which gravity points downwards, 
toward the ground or bottom of the ocean. The vertical coordinate ``z`` 
thus *increases upwards*. We locate the surface at ``z=0``. This means that if
the boundary layer has depth ``h``, the bottom of the boundary layer is 
located at ``z=-h``.

## Governing equations

The one-dimensional, horizontally-averaged boundary-layer equations for 
horizontal momentum ``U`` and ``V``, salinity ``S``, and 
temperature ``T`` are 

```math
\beqs
U_t =   f V - \d_z \overline{w u} - F^u \c \label{xmomentum} \\
V_t = - f U - \d_z \overline{w v} - F^v \c \\
S_t =       - \d_z \overline{w s} - F^S \c \\
T_t =       - \d_z \overline{w \theta} - F^T \c \label{temperature}
\eeqs
```

where subscripts ``t`` and ``z`` denote derivatives with respect to time 
and the vertical coordinate ``z`` and ``f`` is the Coriolis parameter. 
The lowercase variables ``u``, ``v``, ``s``, and ``\theta`` refer to the 
three-dimensional perturbations from horizontal velocity, salinity, and 
temperature, respectively. 
In \eqref{xmomentum}--\eqref{temperature}, internal and boundary forcing of a
quantity ``\Phi`` is denoted ``F^\Phi``, while vertical turbulent fluxes are

```math
G^\Phi = \overline{w \phi} \p
```

# Numerics

One of the surprising aspects of boundary layer parameterization is that 
spatial discretization and time-stepping considerations are not less important
or even distinct from physical aspects of approximation and modeling.

## Spatial discretization

`OceanBoundaryLayerModels.jl` uses a finite volume method to discretize the
oceanic boundary layer in ``z``.

## Time-stepping

Models for the oceanic boundary layer are partial differential equations of 
the form

```math
\beq
\Phi_t = L \Phi + N(\Phi) \c
\label{mathematicalform}
\eeq
```

where ``\Phi`` is a variable like velocity, temperature, or salinity, ``L`` is 
a linear operator, and ``N`` is a nonlinear operator. The linear and nonlinear
parts of the boundary layer model are delineated to permit hybrid 
implicit/explicit time-stepping schemes.

### Forward-Euler method

An explicit forward Euler time integration scheme discretizes
\eqref{mathematicalform} in time with

```math
\beq
\Phi^{n+1} = \Phi^{n} + \Delta t \left [ L \Phi^n + N(\Phi^n) \right ] \,
\eeq 
```

where the superscripts ``n`` and ``n+1`` denote the solution at 
time-step ``n`` and ``n+1``, respectively.
