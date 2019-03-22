# Physics of the oceanic boundary layer

```math
\newcommand{\c}{\, ,}

\newcommand{\r}[1]{\mathrm{#1}}

\newcommand{\ee}{\mathrm{e}}

\newcommand{\beq}{\begin{equation}}
\newcommand{\eeq}{\end{equation}}

\newcommand{\beqs}{\begin{gather}}
\newcommand{\eeqs}{\end{gather}}
```

The dynamics of the ocean's boundary layer are dominated by forcing from the
atmosphere. Atmospheric forcing includes momentum forcing by winds, salinity
forcing by evaporation and precipitation, and heat forcing by latent heat
fluxes, sensible heat fluxes, and incoming and outgoing radiation.

# Models of the oceanic boundary layer

Models for the oceanic boundary layer are partial
differential equations that approximate the effects of

* Atmospheric and radiative forcing;
* Parameterization of convection due to surface cooling;
* Parameterization of mechanical mixing by mixed layer turbulence due mainly
  to wind forcing of boundary-layer currents.

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
u_t - f v = -G^u_z - F^u_z \c \\
v_t + f u = -G^v_z - F^v_z \c \\
      S_t = -G^S_z - F^S_z \c \\
      T_t = -G^T_z - F^T_z \c
\eeqs
```

where subscripts ``t`` and ``z`` denote derivatives with respect to time
and the vertical coordinate ``z``, ``f`` is the Coriolis parameter,
``G^\phi = \overline{w \phi}`` denotes the turbulent vertical flux of
a quantity ``\phi``, while ``F^\phi`` denotes vertical fluxes due to
forcing.

### Temperature forcing

We write the temperature forcing ``F^T`` as

```math
F^T = F^{\r{lat}} + F^{\r{sens}} + F^{\r{longwave}}
        + F^{\r{shortwave}} \c
```

in terms of the four contributions from latent heating, sensible heating,
outgoing longwave radiation, and incoming shortwave radiation.
The first three contributions are implemented as effective boundary conditions
in the uppermost gridpoints of the model.
Shortwave radiation, on the other hand, heats the interior of the boundary
layer.
We parameterize the effect of interior heating by shortwave radiation by
dividing the shortwave spectrum into infrared (IR) and ultraviolet (UV)
components and introducing a ``z``-dependent absorption function such that

```math
\beq
F^{\r{shortwave}}(z) = F^{\r{shortwave}}_0
    \left ( \alpha_{\r{IR}} \exp \left [ z/d_{\r{IR}} \right ]
          + \alpha_{\r{UV}} \exp \left [ z/d_{\r{UV}} \right ] \right ) \c
\label{shortwaverad}
\eeq
```

where ``F^{\r{shortwave}}_0`` is the incoming shortwave radiation at the
surface, which is provided as an input to the boundary layer model.
In \eqref{shortwaverad},
``d_{\r{IR}}`` and ``d_{\r{UV}}`` are the
penetration scales of infrared and ultraviolet radiation, and
``\alpha_{\r{IR}}`` and ``\alpha_{\r{UV}}``
are the fractions of total shortwave radiation with infrared and
ultraviolet wavelengths, respectively, such that
``\alpha_{\r{IR}} + \alpha_{\r{UV}} = 1``.

### Salinity forcing

Salinity forcing ``F^S`` at the surface is

```math
\beq
F^S(z=0) = S(z=0)(E-P) \c
\label{salinityflux}
\eeq
```
where ``S`` is salinity, ``E`` is evaporation, and ``P`` is precipitation.
Equation \eqref{salinityflux} is implemented as an effective boundary
condition in the uppermost grid points of the model.

### Momentum forcing

Momentum forcing at the surface is

```math
\beq
F^u(z=0) = \frac{\tau^x}{\rho_0} \c
\label{uforcing}
\eeq
```

where ``\tau^x`` is the wind stress in the ``x``-direction with units of
``N/\r{m^2}``, and ``\rho_0`` is a reference density.
A similar form is used for ``y``-momentum.
Equation \eqref{uforcing} is implemented as an effective boundary condition.
