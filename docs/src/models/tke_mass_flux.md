# An eddy-diffusivity mass-flux (EDMF) scheme with prognostic turbulent kinetic energy

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
\renewcommand{\F}[2]    {\Upsilon^{#1}_{#2}}
\renewcommand{\C}[2]    {C^{#1}_{#2}}

\newcommand{\uwind}     {u_\star}
\newcommand{\ubuoy}     {w_\star}

\newcommand{\defn}      {\stackrel{\r{def}}{=}}
```

The EDMF family of schemes parameterizes turbulent convection by introducing a
conditional average that partitions the subgrid boundary layer flow into a
turbulent 'environment' with area ``a_0``,
and non-turbulent updrafts and downdrafts with areas ``a_i`` for ``i>0``.

We use an EDMF scheme with a prognostic equation for turbulent kinetic energy.
At the moment, no plume model is implemented, which means our
EDMF scheme only models mixing due to the presence of turbulent kinetic energy
in stress-driven boundary layers.
A description of mass flux schemes is included below for pedagogical purposes.

## Eddy diffusivities and mass fluxes

In all schemes, the turbulent velocity fluxes are parameterized with an
eddy diffusivity.
For the ``x``-velocity ``U``, for example,

```math
\beq
\overline{w u} = \d_z \left ( K_U \d_z U \right ) \c
\eeq
```
where ``K_U`` is the eddy diffusivity for ``x`` and ``y``-momentum.
The turbulent flux of scalars ``\phi`` such as temperature and salinity
is parameterized by both a turbulent flux and mass transport,

```math
\beq
\overline{u \phi} = \d_z \left ( K \d_z \Phi \right )
  - \d_z \sum_i a_i \tilde W_i \tilde \Phi_i \c
\eeq
```
where ``\tilde \Phi_i`` is the difference between the average of
``\phi`` within domain ``i`` and the total horizontal average ``\Phi``:

```math
\begin{align}
\tilde \Phi_i & \defn \left ( \frac{1}{A_i} \int_{A_i} \r{d} A - \frac{1}{A} \int_A \r{d} A \right ) \phi \c \\
&= \Phi_i - \Phi \c
\end{align}
```
where we have introduced the notation ``\Phi_i`` to denote the average of ``\phi`` within
the environment or updraft area ``A_i``.
The terms ``a_i \tilde W_i \tilde \Phi_i`` account for the vertical transport of ``\phi``
by environment and updraft vertical velocities ``W_i``.

## Turbulent kinetic energy equation

A relatively simple EDMF scheme emerges in the limit of vanishing updrafts
and downdrafts, in which case ``W = W_0 = a_0 = 0``.
In the 1.5-order version of this closure, turbulent diffusivity is modeled
via the prognostic turbulent kinetic energy (TKE) equation

```math
\beq
\d_t e = \d_z \left ( K_e \d_z e \right )
    + K_U \left [ \left ( \d_z U \right )^2 + \left ( \d_z V \right )^2 \right ] 
    + \overline{w b} - \C{D}{} \frac{e^{3/2}}{\ell} \c
  \label{TKE}
\eeq
```
where ``\C{D}{}`` is the TKE dissipation parameter, ``K_U`` is the eddy diffusivity for momentum, 
``K_e`` is the eddy diffusivity for turbulent kinetic energy, and ``\overline{wb}`` is the
buoyancy flux,

```math
\beq
\overline{w b} = - g \left ( \alpha K_C \d_z T + \beta K_C \d_z S \right ) \c
\eeq
```

where ``g`` is gravitational acceleration, ``\alpha`` and ``\beta`` are the thermal 
expansion and haline contraction coefficients, and ``K_C`` is the eddy diffusivity for tracers.
Eddy diffusivities are defined

```math
\beq
K_\Phi = C^K_\Phi \ell_\Phi \sqrt{e}
\eeq
```

where ``C^K_\Phi`` is a model constant and ``\ell_\Phi`` is the mixing length for 
the quantity ``\Phi``.
