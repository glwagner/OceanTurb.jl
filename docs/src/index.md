```math
\newcommand{\c}{\, ,}

\newcommand{\r}[1]{\mathrm{#1}}

\newcommand{\ee}{\mathrm{e}}

\newcommand{\beq}{\begin{equation}}
\newcommand{\eeq}{\end{equation}}

\newcommand{\beqs}{\begin{gather}}
\newcommand{\eeqs}{\end{gather}}
```


# OceanBoundaryLayerModels.jl

`OceanMixedLayers.jl` implements models for the ocean's boundary layer.

# Physics of the oceanic boundary layer

The dynamics of the ocean's boundary layer are dominated by forcing from the
atmosphere. Atmospheric forcing includes momentum forcing by winds, salinity
forcing by evaporation and precipitation, and heat forcing by latent heat
fluxes, sensible heat fluxes, and incoming and outgoing radiation.

## Governing equations for 1D boundary layers

We write the one-dimensional boundary-layer balances for horizontal momentum, 
salinity, and temperature as

```math
\beqs
u_t - f v = -G^u_z - F^u_z \c \\
v_t + f u = -G^v_z - F^v_z \c \\
      S_t = -G^S_z - F^S_z \c \\
      T_t = -G^T_z - F^T_z \c
\eeqs
```

where ``G^\phi = \overline{w \phi}`` denotes the turbulent vertical flux of 
a quantity ``\phi``, while ``F^\phi`` denotes vertical fluxes due to 
forcing.

### Temperature forcing

We write the temperature forcing ``F^T`` as

```math
F^T = F^{\r{lat}} + F^{\r{sens}} + F^{\r{longwave}} 
        + F^{\r{shortwave}} \c
```

in terms of the four contributions from latent heating, sensible heating, 
outgoing longwave radiation, and incoming shortwave radiation. The first three
contributions are implemented as boundary conditions. Shortwave 
radiation, on the other hand, heats the interior of the boundary layer. We 
parameterize the effect of interior heating by shortwave radiation by
dividing the specrtum into infrared (IR) and ultraviolet (UV)
components and introducing an absorption function

```math
\beq
F^{\r{shortwave}} = F^{\r{shortwave}}_0
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
