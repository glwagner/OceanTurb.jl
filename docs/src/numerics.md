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

An ASCII-art respresentation of an example grid with `N=3` is

```text
 ▲ z
 |
         i=4           *         
                j=4   ===  Top   ▲              
         i=3           *         | Δf (i=3)
                j=3   ---        ▼
         i=2           *             ▲            
                j=2   ---            | Δc (j=2)
         i=1           *             ▼  
                j=1   ===  Bottom
         i=0           *           
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
The (negative of the) divergence of the diffusive flux at node ``i`` is therefore

```math
\newcommand{\Kdz}[1]{K_{#1} \left ( \d_z \phi \right )_{#1} }
\begin{align}
\left ( \d_z K \d_z \phi \right )_i &= \frac{ \Kdz{i+1} - \Kdz{i} }{\Delta f_i} \c \\
&= \frac{
          \tfrac{K_{i+1}}{\Delta c_{i+1}} \phi_{i+1}
        - \left ( \tfrac{K_{i+1}}{\Delta c_{i+1}} + \tfrac{K_i}{\Delta c_i} \right ) \phi_i
         + \tfrac{K_i}{\Delta c_i} \phi_{i-1}}{\Delta f_i} \p
\end{align}
```

In the top cell where ``i=N``, the diffusive flux is

```math
\begin{align}
\left ( \d_z K \d_z \phi \right )_N &= \frac{ - F_{\mathrm{top}} - K_N \left ( \d_z \phi \right )_N}{\Delta f_N} \c \\
&= -\frac{F_{\mathrm{top}}}{\Delta f_N} - \frac{K_N \phi_N - K_N \phi_{N-1}}{\Delta f_N \Delta c_N} \p
\end{align}
```

In the bottom cell where ``i=1``, on the other hand, the diffusive flux is

```math
\begin{align}
\left ( \d_z K \d_z \phi \right )_1 &= \frac{  K_2 \left ( \d_z \phi \right )_2 + F_{\mathrm{bottom}}}{\Delta f_1} \c \\
&= \frac{F_{\mathrm{bottom}}}{\Delta f_1} + \frac{K_2 \phi_2 - K_2 \phi_1}{\Delta f_1 \Delta c_2}
\end{align}
```


# Time-stepping

To integrate ocean surface boundary layer models forward in time,
we implement various explicit and implicit-explicit time-stepping schemes.

## Explicit schemes

Our explicit time-stepping schemes integrate partial differential equations
of the form

```math
\beq
\d_t \Phi = R(\Phi) \c
\label{explicitform}
\eeq
```

where ``\Phi`` is a variable like velocity, temperature, or salinity and ``R``
is an arbitrary function representing any number of processes, including turbulent
diffusion and internal forcing.

## Implicit-explicit time-stepping

Our mixed implicit-explicit time-stepping schemes integrate partial differential equations
of the form
```math
\beq \label{implicitdiffusion}
\d_t \Phi - \d_z K \d_z \phi = R(\Phi) \c
\eeq
```
where ``K`` is a diffusivity that,  in general, depends on the model state
and thus the time-step.
These implicit-explicit schemes treat the diffusive term on the left
of \eqref{implicitdiffusion} implicitly.

### Forward Euler method

The explicit forward Euler time integration scheme uses
the temporal discretization \eqref{explicitform}:

```math
\beq
\Phi^{n+1} = \Phi^{n} + \Delta t \, R \left ( \Phi^n \right )
\eeq
```

where the superscripts ``n`` and ``n+1`` denote the solution at
time-step ``n`` and ``n+1``, respectively.

### Backward Euler method

The backward Euler method uses the temporal discretization

```math
\beq
\Phi^{n+1} - \Delta t \left ( \d_z K^n \d_z \right ) \Phi^{n+1} = \Phi^n + \Delta t R(\Phi^n)
\eeq
```

The ``z``-derivatives in the diffusive term generate a matrix problem to be solved
for ``\Phi^{n+1}``:

```math
\beq
L_{ij} \Phi^{n+1}_j = \left [ \Phi^n + \Delta t R \left ( \Phi^n \right ) \right ]_i
\eeq
```

where ``L_{ij}`` is a matrix, and the subscripts ``i`` or ``j`` denote grid points
``i`` or ``j``.
Note that the diffusive operator that contributes to ``L_{ij}`` does not include fluxes
across boundary faces; fluxes through boundary faces due either to Dirichlet (Value)
boundary conditions or non-zero fluxes must be included in ``R \left ( \Phi \right )``.
For the diffusive problems considered by our backward Euler solver, ``L_{ij}`` has the form

```math
\beq
L_{ij} = \left [ \begin{matrix}
1 + \Delta t \tfrac{K^n_2}{\Delta f_1 \Delta c_2}
  & -\Delta t \tfrac{K^n_2}{\Delta f_1 \Delta c_2} & \cdot & \cdot & \cdot & \cdot & \cdot \\
- \Delta t \tfrac{K^n_1}{\Delta f_1 \Delta c_1}
  & 1 + \tfrac{\Delta t}{\Delta f_1} \left (\tfrac{K^n_1}{\Delta c_1} + \tfrac{K^n_2}{\Delta f_1 \Delta c_2} \right )
    & -\Delta t \tfrac{K^n_2}{\Delta c_2} & \cdot & \cdot & \cdot & \cdot \\
\cdot & \ddots & \ddots & \ddots & \cdot & \cdot & \cdot \\
\cdot & \cdot
  & - \Delta t \tfrac{K_i}{\Delta c_i \Delta f_i}
  & 1 - \tfrac{\Delta t}{\Delta f_i} \left ( \tfrac{K_{i+1}}{\Delta c_{i+1}} + \tfrac{K_i}{\Delta c_i} \right )
  & - \Delta t \tfrac{K_{i+1}}{\Delta c_{i+1} \Delta f_{i+1}} & \cdot & \cdot \\
\cdot & \cdot & \cdot & \ddots & \ddots & \ddots & \cdot \\
\cdot & \cdot & \cdot & \cdot & \cdot & - \Delta t \tfrac{K^n_N}{\Delta c_N \Delta f_N}
  & 1 + \Delta t \tfrac{K^n_N}{\Delta c_N \Delta f_N}
\end{matrix} \right ]
\eeq
```
