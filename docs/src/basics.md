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

* Internal and surface fluxes of heat, salinity, and momentum due to
  - absorption of incoming solar radiation;
  - cooling by outgoing radiation;
  - latent and sensible heat exchange with the atmosphere;
  - evaporation and precipitation;
  - wind forcing;
* Vertical turbulent fluxes due to
  - convection / gravitational instability;
  - mechanical mixing due to wind and boundary current shear 

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
horizontal momentum, salinity, and temperature are 

```math
\beqs
u_t =   f v - G^u_z - F^u \c \label{xmomentum} \\
v_t = - f u - G^v_z - F^v \c \\
S_t =       - G^S_z - F^S \c \\
T_t =       - G^T_z - F^T \c \label{temperature}
\eeqs
```

where subscripts ``t`` and ``z`` denote derivatives with respect to time 
and the vertical coordinate ``z`` and ``f`` is the Coriolis parameter. 
In \eqref{xmomentum}--\eqref{temperature}, internal and boundary forcing of a
quantity ``\phi`` is denoted ``F^\phi``, while vertical turbulent fluxes are

```math
G^\phi = \overline{w \phi} \p
````

# Numerical modeling of the oceanic boundary layer

## Basic form

Models for the oceanic boundary layer are partial differential equations of 
the form

```math
\beq
\phi_t = L \phi + N(\phi) \c
\label{mathematicalform}
>>>>>>> faf0627972293e4cfa03da9de9e9daa31143d1d5:docs/src/numerics.md
\eeq
```

where ``phi`` is a variable like velocity, temperature, or salinity, ``L`` is 
a linear operator, and ``N`` is a nonlinear operator.

## Time-stepping

An explicit forward Euler time integration scheme discretizes
\eqref{mathematicalform} in time with

```math
\beq
\phi^{n+1} = \phi^{n} + \Delta t \left [ L \phi^n + N(\phi^n) \right ] \,
\eeq 
```

where the superscripts ``n`` and ``n+1`` denote the solution at 
time-step ``n`` and ``n+1``, respectively.
