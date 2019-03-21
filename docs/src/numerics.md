# Numerical methods

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

## Spatial discretization

`OceanTurb.jl` uses a one-dimensional finite-volume method
to discretize momentum, temperature, salinity, and other variables.

An ASCII-art respresentation of an example grid with `nz=3` is

```text
 ▲ z
 |              j=4   ===  Top   ▲              
         i=3           *         | dzf (i=3)
                j=3   ---        ▼
         i=2           *             ▲            
                j=2   ---            | dzc (j=2)
         i=1           *             ▼  
                j=1   ===  Bottom
```

where the double lines indicate the top and bottom of the domain,
the single lines indicate "face" boundaries, the
`i`'s index cell centers (nodes) and `j`'s index the ``z``-location
of cell interfaces (faces).
Horizontal momentum and tracer variables are located at cell centers,
while fluxes of these quantities (and vertical-velocity-like variables when present)
are located at cell faces.

### Derivatives and diffusive fluxes

The derivative of a quantity ``\phi`` at face ``i`` is

```math
\beq
\left( \d_z \phi \right )_i = \frac{\phi_i - \phi_{i-1}}{\Delta c_i} \c
\eeq
```

where ``\phi_i`` denotes the value of ``\phi`` at cell ``i``, and
``\Delta c_i = z_{c, i} - z_{c, i-1}`` is the vertical separation between
node ``i`` and cell point ``i-1``.

With diffusivity given on cell interfaces, the diffusive flux
at face ``i`` is just ``K_i \left ( \d_z \phi \right )_i``.
The divergence of the diffusive flux at node ``i`` is therefore

```math
\newcommand{\Kdz}[1]{K_{#1} \left ( \d_z \phi \right )_{#1} }
\begin{align}
\left ( \d_z K \d_z \phi \right )_i &= \frac{ \Kdz{i+1} - \Kdz{i} }{\Delta f_i} \c \\
&= \frac{
          \tfrac{K_{i+1}}{\Delta c_{i+1}} \phi_{i+1}
        - \left ( \tfrac{K_{i+1}}{\Delta c_{i+1}} + \tfrac{K_i}{\Delta c_i} \right ) \phi_i
         + \tfrac{K_i}{\Delta c_i} \phi_{i-1}}{\Delta f_i}
\end{align}
```

# Time-stepping

To integrate ocean surface boundary layer models forward in time,
we implement explicit time-stepping schemes for partial differential
equations of the form

```math
\beq
\d_t \Phi = R(\Phi) \c
\label{explicitform}
\eeq
```

where ``\Phi`` is a variable like velocity, temperature, or salinity and ``R``
is an arbitrary function representing any number of processes, including turbulent
diffusion and internal forcing.

In the future, we will implement implicit-explicit schemes for partial
differential equations of the form
```math
\beq \label{implicitdiffusion}
\d_t \Phi - \d_z K \d_z \phi = R(\Phi) \c
\eeq
```
where ``K`` is a diffusivity.
These implicit-explicit schemes will treat diffusive terms on the left
of \eqref{implicitdiffusion} implicitly.

### Forward-Euler method

An explicit forward Euler time integration scheme uses
the temporal discretization \eqref{explicitform}:

```math
\beq
\Phi^{n+1} = \Phi^{n} + \Delta t R \left ( \Phi^n \right )
\eeq
```

where the superscripts ``n`` and ``n+1`` denote the solution at
time-step ``n`` and ``n+1``, respectively.
