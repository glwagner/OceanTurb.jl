# The eddy-diffusivity mass-flux (EDMF) schemes

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

\newcommand{\defn}      {\stackrel{\r{def}}{=}}
```

The EDMF family of schemes parameterizes turbulent convection by introducing a
conditional average that partitions the subgrid boundary layer flow into a
turbulent 'environment' with area ``a_0``,
and non-turbulent updrafts and downdrafts with areas ``a_i`` for ``i>0``.

We consider two types of EDMF schemes: those with prognostic equations
that model the time-evolution and spatial distribution of turbulent kinetic
energy (TKE), and those that parameterize the effect of turbulent environmental
mixing with a 'K-profile'.

## Turbulent eddy diffusivty and mass flux

In all schemes, the turbulent velocity fluxes are parameterized with an
eddy diffusivity.
For the ``x``-velocity ``U``, for example,

```math
\beq
\overline{w u} = \d_z \left ( K \d_z U \right ) \c
\eeq
```
where ``K`` is the eddy diffusivity.
We consider various models for eddy diffusivity ranging from a model similar to
the K-profile parameterizaton (KPP), and a formulation
in terms of a prognostic, time- and ``z``-dependent turbulent kinetic energy
variable, ``e``.
The turbulent flux of scalars ``\phi`` such as temperature and salinity
is parameterized by both a turbulent flux and mass transport,

```math
\beq
\overline{u \phi} = \d_z \left ( K \d_z \Phi \right )
  - \d_z \sum_i a_i W^*_i \Phi^*_i \c
\eeq
```
where ``\Phi^*_i`` is the difference between the average of
``\phi`` within domain ``i`` and the total horizontal average ``\Phi``:

```math
\begin{align}
\Phi^*_i &\defn \left ( \frac{1}{A_i} \int_{A_i} \r{d} A - \frac{1}{A} \int_A \r{d} A \right ) \phi \c \\
&= \Phi_i - \Phi \c
\end{align}
```
where we have introduced the notation ``\Phi_i`` to denote the average of ``\phi`` within
the environment or updraft area ``A_i``.
The terms ``a_i W^*_i \Phi^*_i`` account for the vertical transport of ``\phi``
by environment and updraft vertical velocities ``W_i``.

## Zero-plume, 1.5-order EDMF scheme

A relatively simple EDMF scheme emerges in the limit of vanishing updrafts
and downdrafts, in which case ``W = W_0 = a_0 = 0``.
In the 1.5-order version of this closure, turbulent diffusivity is modeled
via the prognostic turbulent kinetic energy (TKE) equation

```math
\beq
\d_t e = K \left [ \left ( \d_z U \right )^2 + \left ( \d_z V \right )^2 \right ] + \d_z \left ( K \d_z e \right )
  - K \d_z B - \C{\ep}{} \frac{e^{3/2}}{\F{\ell}{}} \c
  \label{TKE}
\eeq
```
where ``\C{\ep}{} = 2.0`` is a model parameter,
``\d_z B = g \left( \alpha \d_z T - \beta \d_z S \right )`` is the buoyancy gradient in terms of gravitational acceleration ``g`` and thermal expansion and haline contraction coefficients ``\alpha`` and ``\beta``,
and ``K`` is the eddy diffusivity defined in terms of turbulent 'velocity' ``\sqrt{e}`` and a mixing length ``\F{\ell}{}``:

```math
\beq
K = \C{K}{} \underbrace{
      \C{\kappa}{} z \left ( 1 - \F{a}{} \tfrac{F_b}{\uwind^3} z \right )^{\F{n}{}}}
        _{\defn \F{\ell}{}}
        \, \sqrt{e} \p
        \label{eddydiffusivity}
\eeq
```
In \eqref{eddydiffusivity}, ``F_b = g \left ( \alpha F_\theta - \beta F_s \right )`` is the buoyancy flux define in terms of temperature and salinity fluxes ``F_\theta`` and ``F_s``, and ``\uwind \defn | \b{F}_u |^{1/2}`` is the friction velocity defined in terms of velocity flux ``\b{F}_u = \b{\tau} / \rho_0`` or wind-stress ``\b{\tau}`` and reference density ``\rho_0``.
``\C{\kappa}{} = 0.41`` and ``\C{K}{} = 0.1`` in \eqref{eddydiffusivity} are the 'Von Karman' and eddy diffusivity model parameters, respectively.
``\F{a}{}`` and ``\F{n}{}`` in \eqref{eddydiffusivity} are piecewise constant model functions
that model the effect of boundary layer stability on the mixing length and are

```math
\beq
\F{a}{} = \left \{ \begin{matrix}
-100 & \text{for unstable boundary layers with } F_b > 0 \\
2.7 & \text{for stable boundary layers with } F_b \le 0
\end{matrix} \right . \c
\eeq
```

and

```math
\beq
\F{n}{} = \left \{ \begin{matrix}
0.2 & \text{for unstable boundary layers with } F_b > 0 \\
-1 & \text{for stable boundary layers with } F_b \le 0
\end{matrix} \right . \p
\eeq
```

Note that parameterized buoyancy flux ``\overline{w b} \defn -K \d_z B`` appears in the
TKE equation \eqref{TKE}.
