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
\newcommand{\K}         {\mathrm{KE}}        

\newcommand{\btau}      {\b{\tau}} % wind stress vector

% Model functions and constants
\renewcommand{\F}[2]      {\Upsilon^{#1}_{#2}}
\renewcommand{\C}[2]      {C^{#1}_{#2}}

\newcommand{\uwind}     {\varpi_{\tau}}
\newcommand{\ubuoy}     {\varpi_b}
```

The K-Profile-Parameterization, or "KPP", is proposed by
[Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
as a model for convection- and wind-driven mixing in the upper ocean.
In KPP, vertical turbulent fluxes of a quantity ``\phi`` are parameterized as

```math
\beq
\overline{w \phi} = - K_\Phi \d_z \Phi + N_\phi \c
\eeq
```

where ``\Phi`` is the resolved or horizontally-averaged quantity, ``K_\Phi`` is a
turbulent diffusivity, and ``N_\phi`` is a 'non-local' flux.

The non-local flux and turbulent diffusivity are defined to vanish at the surface, and at the bottom of the 'mixing layer' ``h``, which roughly
corresponds to the depth at which turbulent fluxes and turbulent kinetic energy
decay to zero.

## Description of the model

Below, we denote 'model parameters' as ``\C{\mathrm{label}}{\mathrm{var}}``, where
'label' describes the parameter and 'var' is a variable like ``U, V, T`` or ``S``.

Buoyancy is then defined in terms of ``T`` and ``S`` as

```math
\begin{align}
B & \equiv - \frac{g \rho'}{\rho_0} \\
  &     = g \left [ \alpha \left ( T - T_0 \right ) - \beta \left ( S - S_0 \right ) \right ] \c
\end{align}
```

where ``g = 9.81 \, \mathrm{m \, s^{-2}}, \alpha = 2 \times 10^{-4} \, \mathrm{K^{-1}}``, and ``\beta = 8 \times 10^{-5}``,
are gravitational acceleration, the thermal expansion coefficient, and the
haline contraction coefficient, respectively.
Buoyancy forcing is equivalently

```math
\beq
F_b = g \left ( \alpha F_\theta - \beta F_s \right ) \p
\eeq
```

The turbulent velocity scales associated with buoyancy and wind forcing are

```math
\beq
\omega_b \equiv | h F_b |^{1/3} \qquad \text{and} \qquad \omega_\tau \equiv | \b{\tau} / \rho_0 |^{1/2} \p
\eeq
```

where ``h`` is the depth of the 'mixing layer', or the depth to which
turbulent mixing and turbulent fluxes penetrate, ``\b{\tau}`` is wind stress,
and ``\rho_0 = 1028 \, \mathrm{kg \, m^{-3}}`` is a reference density.

We also define the ratios

```math
\beq
r_b \equiv \left ( \frac{\omega_b}{\omega_\tau} \right )^3 \qquad \text{and}
\qquad r_\tau \equiv \left ( \frac{\omega_\tau}{\omega_b} \right )^3 = \frac{1}{r_b} \p
\eeq
```

### The boundary layer depth, ``h``

The boundary layer depth ``h`` is defined implicitly via the bulk Richardson number criterion

```math
\beq \label{bulk_ri}
\C{\Ri}{} = \frac{h \Delta B(-h)}{| \Delta \b{U}(-h)|^2 + \F{\K}{}(-h)} \p
\eeq
```

where the critical ``\Ri`` is ``\C{\Ri}{} = 0.3``. The operator ``\Delta`` is defined

```math
\beq
\Delta \Phi(z) = -\frac{1}{\C{\ep}{} z} \int_{\C{\ep}{} z}^0 \Phi(z') \di z' - \Phi(z) \c
\eeq
```

where ``\C{\ep}{} = 0.1`` is the surface layer fraction.
The function ``\F{\K}{}(z)`` is

```math
\beq  \label{unresolved_ke}
\F{\K}{}(z) = \C{\K}{} (-z)^{4/3} \sqrt{ \max \left [ 0, B_z(z) \right ] } \max \left [ 0, F_b^{1/3} \right ] \p
\eeq
```

The unresolved kinetic energy constant is ``\C{\K}{} = 4.32``.


where ``g`` is gravitational acceleration and ``\rho_0, T_0, S_0`` are reference densities, temperatures, and salinities.
``\rho'`` is the density deviation from the reference.

### Non-local flux

The non-local flux is defined only for ``T`` and ``S``, and is

```math
\beq
N_\Phi = \C{N}{} F_\Phi d (1 - d)^2 \c
\eeq
```
where ``d = -z/h`` is a non-dimensional depth coordinate and ``\C{N}{} = 6.33``.

### Turbulent Diffusivity

The KPP diffusivity is defined

```math
\beq
K_\Phi = h \F{w}{\Phi} d ( 1 - d )^2 \c
\eeq
```
where ``\F{w}{\Phi}`` is the turbulent velocity scale.
In wind-driven turbulence under stable buoyancy forcing such that ``F_b < 0``, the turbulent velocity scale is

```math
\beq
\F{w}{\Phi} = \frac{ \C{\kappa}{} \uwind}{1 + \C{\mathrm{stab}}{} r_b d} \p
\eeq
```

wherea ``\C{\kappa}{} = 0.4`` is the Von Karman constant and ``\C{\mathrm{stab}}{} = 2.0``.

In wind-driven turbulence but under destabilizing buoyancy forcing, when ``\min \left [ \C{\ep}{}, d \right ] < \C{d}{\phi} r_\tau``,
the turbulent velocity scale is

```math
\beq
\F{w}{\Phi} = \C{\kappa}{} \omega_\tau \left ( 1 + \C{\mathrm{unst}}{} r_b \min \left [ \C{\ep}{}, d \right ] \right )^{n_\Phi} \c
\eeq
```

where ``n_U = 1/4`` for velocities, ``n_T = 1/2`` for tracers, ``\C{\mathrm{unst}}{} = 6.4``, the
transition parameter for velocities is ``\C{d}{U} = 0.5``, and the transition parameter
for tracers is ``\C{d}{T} = 2.5``.

In convection-driven turbulence affected by wind mixing, when ``\min \left [ \C{\ep}{}, d \right ] >= \C{d}{\phi} r_\tau``,
the turbulent velocity scale is

```math
\beq
\F{w}{\Phi} = \C{b}{\Phi} \omega_b \left ( \min \left [ \C{\ep}{}, d \right ] + \C{\tau}{\Phi} r_\tau \right )^{1/3} \c
\eeq
```

where ``\C{b}{U} = \C{b}{V} = 0.215``, ``\C{b}{T} = \C{b}{S} = 2.53``, ``\C{\tau}{U} = \C{\tau}{V} = 0.0806``, and
``\C{\tau}{T} = \C{\tau}{S} = 1.85``.


## Selected tests

See `/test/runtests.jl` for more tests.

### Linear temperature profile and no velocity field

Some simple tests can be defined when the model state is ``U=V=S=0`` and ``T = \gamma z``.
In this case the buoyancy becomes ``B = g \alpha \gamma z`` and the buoyancy gradient is ``B_z = g \alpha \gamma``.
If we further take ``\C{\ep}{} \to 0`` and ``F_b > 0``, and note that the value of ``T``
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
\F{\K}{}(-h) = \C{\K}{} h^{4/3} \sqrt{g \alpha \gamma} F_b^{1/3} \p
\eeq
```

The bulk Richardson number criterion then becomes

```math
\begin{align}
\C{\Ri}{} &= \frac{h \Delta B(-h)}{\F{\K}{}(-h)} \c \\
          &= \frac{g \alpha \gamma h^2 - h B_N}{\C{\K}{} h^{4/3} \sqrt{g \alpha \gamma} F_b^{1/3}} \c \\
          &= \frac{g \alpha \gamma h - B_N}{\C{\K}{} \sqrt{g \alpha \gamma} \left ( h F_b \right )^{1/3}} \c \\
\end{align}
```

Modifying the temperature profile so that ``T_N = B_N = 0`` allows us to
analytically calculate the mixed layer depth:

```math
\beq
h = \left ( \C{\Ri}{} \C{\K}{} \right )^{3/2} \sqrt{F_b} \left ( g \alpha \gamma \right )^{-3/4} \p
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
Setting ``h = 9``, ``\C{\Ri}{}=1``, ``T_0 = 1``, and ``U_0=3`` yields a consistent solution.

### Limiting cases for turbulent velocity scales

Under zero momentum forcing, the turbulent vertical velocity scale is

```math
\beq
\F{w}{\Phi} = \C{b}{\Phi} \left ( \C{\ep}{} \right )^{1/3} | h F_b |^{1/3} \p
\label{buoyscaletest}
\eeq
```

we write the test in \eqref{buoyscaletest} using the depth in \eqref{analyticaldepth}

Under zero buoyancy forcing, the turbulent velocity scale is

```math
\beq
\F{w}{\Phi} = \C{\kappa}{} \omega_\tau \p
\label{windscaletest}
\eeq
```


# Table of model parameters

The model parameters in KPP are

|   Parameter   | Value | Reference              |
|   :-------:   | :---: | ---------              |
| ``\C{\ep}{}`` | 0.1   | pretty much everywhere |


# References

* [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
* [CVMix documentation](https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf)
* [Van Roekel et al (2018)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001336)
