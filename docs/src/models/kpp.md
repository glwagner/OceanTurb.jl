# The K-Profile-Parameterization (KPP)

```math
\newcommand{\c}         {\, ,}
\newcommand{\p}         {\, .}
\newcommand{\d}         {\partial}
\newcommand{\r}[1]      {\mathrm{#1}}
\newcommand{\b}[1]      {\boldsymbol{#1}}
\newcommand{\ee}        {\mathrm{e}}
\newcommand{\di}        {\, \mathrm{d}}
\newcommand{\ep}        {\epsilon}

\newcommand{\beq}       {\begin{equation}}
\newcommand{\eeq}       {\end{equation}}
\newcommand{\beqs}      {\begin{gather}}
\newcommand{\eeqs}      {\end{gather}}

% Non-dimensional numbers
\newcommand{\Ri}        {\mathrm{Ri}}
\newcommand{\SL}        {\mathrm{SL}}
\newcommand{\K}         {\mathcal{E}}
\newcommand{\W}         {\mathcal{W}}

\newcommand{\btau}      {\b{\tau}} % wind stress vector

% Model functions and constants
\renewcommand{\F}[2]      {\Upsilon^{#1}_{#2}}
\renewcommand{\C}[2]      {C^{#1}_{#2}}

\newcommand{\uwind}     {u_\star}
\newcommand{\ubuoy}     {w_\star}

\newcommand{\NL}        {NL}
```

The K-Profile-Parameterization, or "KPP", is proposed by
[Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
as a model for convective and wind-driven mixing in the ocean surface boundary layer.
In KPP, vertical turbulent fluxes of a quantity ``\phi`` are parameterized as

```math
\beq
\overline{w \phi} = - K_\phi \d_z \Phi + \NL_\phi \c
\eeq
```

where ``\Phi`` is the resolved or horizontally-averaged quantity, ``K_\phi`` is a
turbulent diffusivity, and ``NL_\phi`` is a 'non-local' flux.

The non-local flux and turbulent diffusivity are defined to vanish at the surface, and at the bottom of the 'mixing layer' ``h``, which roughly
corresponds to the depth at which turbulent fluxes and turbulent kinetic energy
decay to zero.

The various numerical implementations of 'KPP' have effectively resulted in
the proliferation of practically distinct KPP models in different codes.
Here we describe the implementation of KPP in `OceanTurb`, which is meant to
follow the version of KPP described by the
[Community Vertical Mixing Project (CVMix)](https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf).

The KPP model has three distinct parts:

1. A model for the mixing layer depth, ``h``.
2. A model for the non-local flux, ``\NL_\Phi``.
3. A model for the local diffusivity, ``K_\Phi``.

Below, we denote 'model parameters' as ``\C{\mathrm{label}}{\mathrm{var}}``, where
'label' describes the parameter and 'var' is a variable like ``U, V, T`` or ``S``.

## Mixing depth model in CVMix KPP

The mixing layer depth ``h`` is defined implicitly via the bulk Richardson number criterion

```math
\beq \label{bulk_ri}
\C{\Ri}{} = \frac{h \left ( 1 - \tfrac{1}{2} \C{\ep}{} \right ) \Delta B(-h)}{| \Delta \b{U}(-h)|^2 + \K(-h)} \
\eeq
```

where the critical ``\Ri`` is ``\C{\Ri}{} = 0.3``. The operator ``\Delta`` is defined

```math
\beq
\Delta \Phi(z) = -\frac{1}{\C{\SL}{} z} \int_{\C{\SL}{} z}^0 \Phi(z') \di z' - \Phi(z) \c
\eeq
```

where ``\C{\SL}{} = 0.1`` is the surface layer fraction.
The function ``\K(-h)`` parameterizes unresolved kinetic energy associated with convective
plumes,

```math
\beq  \label{unresolved_ke}
\K(-h) = \C{\K}{} (-z)^{4/3} \sqrt{ \max \left [ 0, B_z(z) \right ] } \max \left [ 0, Q_b \right ]^{1/3} + \C{\K_0}{} \p
\eeq
```

The unresolved kinetic energy constant is ``\C{\K}{} = 4.32`` and the minimum unresolved kinetic energy is ``\C{\K_0}{} = 10^{-11}``.
To solve \eqref{bulk\_ri} for ``h``, we evaluate the right side of \eqref{bulk\_ri} for ``z < 0`` at increasing depths until the right side rises above the critical ``\C{\Ri}{}``.
We then linearly interpolate to find ``h``.

## 'Countergradient' non-local flux model in CVMix KPP

The non-local flux is defined only for ``T`` and ``S``, and is

```math
\beq
NL_\phi = \C{\NL}{} Q_\phi d (1 - d)^2 \c
\eeq
```
where ``d = -z/h`` is a non-dimensional depth coordinate and ``\C{\NL}{} = 6.33``.

## ``K``-Profile model in CVMix KPP

The KPP diffusivity is

```math
\beq
K_\phi = h \W_\Phi(d) d ( 1 - d )^2 \c
\eeq
```

where ``\W_\Phi(d)`` is the turbulent velocity scale for variable ``\phi``.

We define the ratios

```math
\beq
r_b \equiv \left ( \frac{\ubuoy}{\uwind} \right )^3 \qquad \text{and}
\qquad r_\tau \equiv \left ( \frac{\uwind}{\ubuoy} \right )^3 = \frac{1}{r_b} \p
\eeq
```

In wind-driven turbulence under stable buoyancy forcing such that ``Q_b < 0``, the turbulent velocity scale is

```math
\beq
\W_\Phi = \frac{ \C{\tau}{} \uwind}{ \left ( 1 + \C{\mathrm{stab}}{} r_b d \right )^{\C{n}{}}} \p
\eeq
```

where ``\C{\tau}{} = 0.4`` is the Von Karman constant, ``\C{\mathrm{stab}}{} = 2.0``, and ``\C{n}{}=1``.

In wind-driven turbulence but under destabilizing buoyancy forcing, when ``\min \left [ \C{\SL}{}, d \right ] < \C{d}{\phi} r_\tau``,
the turbulent velocity scale is

```math
\beq
\W_\Phi(d) = \C{\tau}{} \uwind \left ( 1 + \C{\mathrm{unst}}{} r_b \min \left [ \C{\SL}{}, d \right ] \right )^{\C{m\tau}{\Phi}} \c
\eeq
```

where ``\C{m\tau}{U}= 1/4``, ``\C{m\tau}{T} = 1/2``, ``\C{\mathrm{unst}}{} = 6.4``, the
transition parameter for velocities is ``\C{d}{U} = 0.5``, and the transition parameter
for tracers is ``\C{d}{T} = 2.5``.

In convection-driven turbulence affected by wind mixing, when ``\min \left [ \C{\SL}{}, d \right ] >= \C{d}{\phi} r_\tau``,
the turbulent velocity scale is

```math
\beq
\W_\Phi(d) = \C{b}{\Phi} \ubuoy \left ( \min \left [ \C{\SL}{}, d \right ] + \C{\tau b}{\Phi} r_\tau \right )^{\C{mb}{\Phi}} \c
\eeq
```

where
``\C{b}{U} = \C{b}{V} = 0.215``,
``\C{b}{T} = \C{b}{S} = 2.53``,
``\C{mb}{U} = \C{mb}{T}=1/3``,
``\C{\tau b}{U} = \C{\tau b}{V} = 0.374``, and
``\C{\tau b}{T} = \C{\tau b}{S} = -0.717``.


## Selected tests

See `/test/test_kpp.jl` for more tests.

### Linear temperature profile and no velocity field

Some simple tests can be defined when the model state is ``U=V=S=0`` and ``T = \gamma z``.
In this case the buoyancy becomes ``B = g \alpha \gamma z`` and the buoyancy gradient is ``B_z = g \alpha \gamma``.
If we further take ``\C{\SL}{} \to 0`` and ``Q_b > 0``, and note that the value of ``T``
in the top grid cell is ``T_N = -\gamma \Delta z / 2``, where ``N`` is the number of
grid points, we find that

```math
\beq
\Delta T(-h) = \gamma h - T_N \c
\eeq
```

and

```math
\beq
\Delta B(-h) = g \alpha \gamma h - B_N \p
\eeq
```

The unresolved kinetic energy function is

```math
\beq
\K(-h) = \C{\K}{} h^{4/3} \sqrt{g \alpha \gamma} Q_b^{1/3} \p
\eeq
```

The bulk Richardson number criterion then becomes

```math
\begin{align}
\C{\Ri}{} &= \frac{h \Delta B(-h)}{\K(-h)} \c \\
          &= \frac{g \alpha \gamma h^2 - h B_N}{\C{\K}{} h^{4/3} \sqrt{g \alpha \gamma} Q_b^{1/3}} \c \\
          &= \frac{g \alpha \gamma h - B_N}{\C{\K}{} \sqrt{g \alpha \gamma} \left ( h Q_b \right )^{1/3}} \c \\
\end{align}
```

Modifying the temperature profile so that ``T_N = B_N = 0`` allows us to
analytically calculate the mixed layer depth:

```math
\beq
h = \left ( \C{\Ri}{} \C{\K}{} \right )^{3/2} \sqrt{Q_b} \left ( g \alpha \gamma \right )^{-3/4} \p
\label{analyticaldepth}
\eeq
```

### Linear temperature profile and linear shear

Consider a piecewise-constant temperature profile

```math
\beq
T(z) =  \left \{ \begin{matrix}
  T_0 & \quad \text{for } z > -h \c \\
  -T_0 & \quad \text{for } z < -h \c
  \end{matrix} \right .
\eeq
```

and a velocity profile

```math
\beq
U(z) =  \left \{ \begin{matrix}
  U_0 & \quad \text{for } z > -h \c \\
  -U_0 & \quad \text{for } z < -h \c
  \end{matrix} \right .
\eeq
```

so that ``T(z=-h) = U(z=-h) = 0``.
We then have ``\Delta T(-h) = T_0`` and ``\Delta U^2(-h) = U_0^2``, so that
with ``g=\alpha=1``,

```math
\beq
\C{\Ri}{} = \frac{h T_0}{U_0^2} \p
\label{sheardepth}
\eeq
```

setting ``h = 9``, ``\C{\Ri}{}=1``, ``T_0 = 1``, and ``U_0=3`` yields a consistent solution.

### Limiting cases for turbulent velocity scales

Under zero momentum forcing, the turbulent vertical velocity scale is

```math
\beq
\W_\phi = \C{b}{\phi} \left ( \C{\ep}{} \right )^{1/3} | h Q_b |^{1/3} \p
\label{buoyscaletest}
\eeq
```
We write the test in \eqref{buoyscaletest} using the depth in \eqref{analyticaldepth}

Under zero buoyancy forcing, the turbulent velocity scale is

```math
\beq
\W_\phi = \C{\tau}{} \uwind \p
\label{windscaletest}
\eeq
```

# Table of model parameters

The default values for adjustable model parameters in KPP are

|   Parameter             | Value       | Description |
|   :-------:             | :---:       | :---:       |
| ``\C{\Ri}{}``           | 0.3         | Bulk Richardson number criterion |
| ``\C{\SL}{}``           | 0.1         | Surface layer fraction |
| ``\C{\K}{}``            | 3.19        | Unresolved kinetic energy constant |
| ``\C{\NL}{}``           | 6.33        | Non-local flux proportionality constant |
| ``\C{\tau}{}``          | 0.4         | Wind mixing constant / von Karman parameter |
| ``\C{\mathrm{stab}}{}`` | 2.0         | Proportionality constant for effect of stable buoyancy forcing on wind mixing |
| ``\C{n}{}``             | 1.0         | Exponent for effect of stable buoyancy forcing on wind mixing |
| ``\C{\mathrm{unst}}{}`` | 6.4         | Proportionality constant for effect of unstable buoyancy forcing on wind mixing |
| ``\C{m\tau}{U}``        | 0.25        | Exponent for effect of unstable buoyancy forcing on wind mixing of momentum |
| ``\C{m\tau}{T}``        | 0.5         | Exponent for effect of unstable buoyancy forcing on wind mixing of momentum |
| ``\C{b}{U}``            | 0.599       | Convective mixing constant for momentum |
| ``\C{b}{T}``            | 1.36        | Convective mixing constant for scalars |
| ``\C{d}{U}``            | 0.5         | Transitional normalized depth for unstable mixing of momentum |
| ``\C{d}{T}``            | 2.5         | Transitional normalized depth for unstable mixing of scalars |
| ``\C{mb}{U}``           | 0.33        | Exponent for effect of wind on convective mixing of momentum |
| ``\C{mb}{T}``           | 0.33        | Exponent for effect of wind on convective mixing of scalars |
| ``K_{u0}``              | ``10^{-5}`` | Interior/background turbulent diffusivity for momentum |
| ``K_{T0}``              | ``10^{-5}`` | Interior/background turbulent diffusivity for temperature |
| ``K_{S0}``              | ``10^{-5}`` | Interior/background turbulent diffusivity for salinity |

Note: all parameters are greater than 0, and ``0 \ge \C{\SL}{} \ge 1``.

The default values for 'non-adjustable' parameters in KPP are

|   Parameter       | Value  | Description |
|   :-------:       | :---:  | :-: |
| ``\C{\tau b}{U}`` | 0.374  | | 
| ``\C{\tau b}{T}`` | -0.717 | |
| ``\C{\K_0}{}``    | 1e-11  | |


# References

* [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
* [CVMix documentation](https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf)
* [Van Roekel et al (2018)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001336)
