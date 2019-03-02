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

`OceanTurb.jl` uses 
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

`OceanTurb.jl` uses a finite-volume method with one-dimensional analog of 
the staggered Arakawa C-grid to discretize momentum, temperature, salinity,
and other variables.

An ASCII-art respresentation of an exmaple grid with `nz=3` is

```text
      ▲ z 
      |   
        
                j=4   ===       ▲              
         i=3           *        | dzf (i=3)
                j=3   ---       ▼
         i=2           *    ▲            
                j=2   ---   | dzc (j=2) 
         i=1           *    ▼  
                j=1   ===     
```

where the double lines indicate the top and bottom of the domain,
the single lines indicate "face" boundaries, the
`i`'s index cell centers and `j`'s index face-located grid points 
and variables.
Horizontal momentum and tracer variables are located at cell centers, 
while fluxes of these quantities (and vertical-velocity-like variables when present)
are located at cell faces. 


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
