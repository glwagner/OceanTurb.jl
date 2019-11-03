# Numerical methods in `OceanTurb.jl`

`OceanTurb.jl` uses a one-dimensional finite-volume method
to discretize momentum, temperature, salinity, and other variables
in the ``z``-direction.
A variety of explicit and implicit-explicit schemes are
implemented for temporal integration.

# Spatial discretization

An ASCII-art respresentation of an example grid with `N=3` is

```text
 ▲ z
 |
         i=4           *         
                j=4   ===  Top   ▲              
         i=3           *         | Δf[3]
                j=3   ---        ▼
         i=2           *             ▲            
                j=2   ---            | Δc[2]
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
The cells at ``i=0`` and ``i=4`` are 'ghost cells', whose values are set
according to the boundary condition.
For a no flux or zero gradient boundary condition, for example, we
would set `c[0]=c[1]` and `c[4]=c[3]`.

### Finite volume derivatives and fluxes

The derivative of a quantity ``\Phi`` at face ``i`` is

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

\beq
\left( \d_z \Phi \right )_i = \frac{\Phi_i - \Phi_{i-1}}{\Delta c_i} \c
\eeq
```

where ``\Phi_i`` denotes the value of ``\Phi`` at cell ``i``, and
``\Delta c_i = z_{c, i} - z_{c, i-1}`` is the vertical separation between
node ``i`` and cell point ``i-1``.

With diffusivities defined on cell interfaces, the diffusive flux
across face ``i`` is ``K_i \left ( \d_z \Phi \right )_i``.
The (negative of the) divergence of the diffusive flux at node ``i`` is therefore

```math
\newcommand{\Kdz}[1]{K_{#1} \left ( \d_z \Phi \right )_{#1} }
\begin{align}
\left ( \d_z K \d_z \Phi \right )_i &= \frac{ \Kdz{i+1} - \Kdz{i} }{\Delta f_i} \c \\
&= \frac{
          \tfrac{K_{i+1}}{\Delta c_{i+1}} \Phi_{i+1}
        - \left ( \tfrac{K_{i+1}}{\Delta c_{i+1}} + \tfrac{K_i}{\Delta c_i} \right ) \Phi_i
         + \tfrac{K_i}{\Delta c_i} \Phi_{i-1}}{\Delta f_i} \p
 \label{fluxdivop}
\end{align}
```

In the top cell where ``i=N``, the diffusive flux is

```math
\begin{align}
\left ( \d_z K \d_z \Phi \right )_N &= \frac{ - Q_{\mathrm{top}} - K_N \left ( \d_z \Phi \right )_N}{\Delta f_N} \c \\
&= -\frac{Q_{\mathrm{top}}}{\Delta f_N} - \frac{K_N \Phi_N - K_N \Phi_{N-1}}{\Delta f_N \Delta c_N} \p
 \label{fluxdivop_top}
\end{align}
```

In the bottom cell where ``i=1``, on the other hand, the diffusive flux is

```math
\begin{align}
\left ( \d_z K \d_z \Phi \right )_1 &= \frac{  K_2 \left ( \d_z \Phi \right )_2 + Q_{\mathrm{bottom}}}{\Delta f_1} \c \\
&= \frac{Q_{\mathrm{bottom}}}{\Delta f_1} + \frac{K_2 \Phi_2 - K_2 \Phi_1}{\Delta f_1 \Delta c_2}
\label{fluxdivop_bottom}
\end{align}
```

For negative advective mass fluxes ``M`` defined at cell centers (corresponding
to downdrafts or down-travelling plumes) imply an advective flux divergence

```math
\beq
\d_z \left ( M \Phi \right )_i = \frac{M_{i+1} \Phi_{i+1} - M_i \Phi_i}{\Delta c_{i+1}}
\eeq
```

# Time integration

To integrate ocean surface boundary layer models forward in time,
we implement various explicit and implicit-explicit time-stepping schemes.
The function `iterate!(model, Δt, Nt)` steps a model forward in time.

Timesteppers in `OceanTurb.jl` integrate equations of the form

```math
\beq \label{equationform}
\d_t \Phi = - \d_z \left ( M \Phi \right ) + \left ( \d_z K \d_z \right ) \Phi - L \Phi + R(\Phi) \c
\eeq
```
where ``\Phi(z, t)`` is a variable like velocity, temperature, or salinity,
``M``, ``K``, and ``L`` are an advective 'mass flux', diffusivity,
and damping coefficient
which are a general nonlinear functions of ``\Phi``, ``z``, and external parameters,
and ``R`` is an arbitrary function representing any number of processes,
including the Coriolis force or external forcing.

## Time integration methods

We implement `iterate!` functions and types for:

* explicit forward Euler
* semi-implicit backward Euler

### Forward Euler method

The explicit forward Euler time integration scheme
obtains ``\Phi`` at time-step ``n+1`` using the formula

```math
\beq
\Phi^{n+1} = \Phi^{n} + \Delta t \, \big [ - \d_z \left ( M^n \Phi^n \right ) + \left ( \d_z K^n \d_z \right ) \Phi^n - L^n \Phi^n + R \left ( \Phi^n \right ) \big ]
\eeq
```

where the superscripts ``n`` and ``n+1`` denote the solution at
time-step ``n`` and ``n+1``, respectively.

### Backward Euler method

The backward Euler method
obtains ``\Phi`` at time-step ``n+1`` using the formula

```math
\beq
\Phi^{n+1}
  + \Delta t \left [ \d_z \left ( M^n \Phi^{n+1} \right ) - \left ( \d_z K^n \d_z \right ) \Phi^{n+1} + L^n \Phi^{n+1} \right ]
    = \Phi^n + \Delta t R \left ( \Phi^n \right )
\eeq
```

The ``z``-derivatives in the advective and diffusive terms generate
an elliptic problem to be solved for ``\Phi^{n+1}`` at each time-step.
In the finite volume discretization used by `OceanTurb.jl`, this elliptic
problem becomes a matrix problem of the form

```math
\beq
\mathcal{L}^n_{ij} \Phi^{n+1}_j = \left [ \Phi^n + \Delta t R \left ( \Phi^n \right ) \right ]_i
\eeq
```

where ``\mathcal{L}^n_{ij}`` is a matrix operator at time-step ``n``, and the subscripts ``i`` or ``j`` denote grid points
``i`` or ``j``.
For the diffusive problems considered by our backward Euler solver, the matrix
multiplication ``\mathcal{L}^n_{ij} \Phi_j^{n+1}`` has the form


```math
\begin{align}
\mathcal{L}^n_{ij} \Phi_j^{n+1} &= \left [ 1 + \Delta t \left ( \d_z M^n - \d_z K^n \d_z + L^n \right ) \right ]_{ij} \Phi_j^{n+1} \\

&= \left [ \begin{matrix}

1 + \Delta t \left ( L^n_1 + \tfrac{K^n_2}{\Delta f_1 \Delta c_2} - \tfrac{M_1^n}{\Delta c_2} \right )
  & \Delta t \left ( \tfrac{M_2^n}{\Delta c_2} - \tfrac{K^n_2}{\Delta f_1 \Delta c_2} \right )
    & \cdot & \cdot & \cdot & \cdot \\

\ddots & \ddots & \ddots & \cdot & \cdot & \cdot \\

\cdot
  & - \Delta t \tfrac{K^n_i}{\Delta c_i \Delta f_i}
  & 1 + \Delta t
    \left [ L^n_i + \tfrac{K^n_{i+1}}{\Delta f_i \Delta c_{i+1}} + \tfrac{K^n_i}{\Delta f_i \Delta c_i} - \frac{M^n_i}{\Delta c_{i+1}} \right ]
  & \Delta t \left ( \tfrac{M^n_{i+1}}{\Delta c_{i+1}} - \tfrac{K^n_{i+1}}{\Delta c_{i+1} \Delta f_{i+1}} \right ) & \cdot & \cdot \\

\cdot & \cdot & \ddots & \ddots & \ddots & \cdot \\

\cdot & \cdot & \cdot & \cdot
  & - \Delta t \tfrac{K^n_N}{\Delta c_N \Delta f_N}
  & 1 + \Delta t \left ( L^n_N + \tfrac{K^n_N}{\Delta c_N \Delta f_N} - \tfrac{M^n_N}{\Delta c_{N+1}} \right )

\end{matrix} \right ]

\left [ \begin{matrix}
\Phi^{n+1}_1 \\[1.1ex]
\vdots \\[1.1ex]
\Phi^{n+1}_i \\[1.1ex]
\vdots \\[1.1ex]
\Phi^{n+1}_N
\end{matrix} \right ]
\label{implicitoperatormatrix}
\end{align}
```

To form the matrix operator in \eqref{implicitoperatormatrix}, we have used
the second-order flux divergence finite difference operators in
\eqref{fluxdivop}--\eqref{fluxdivop_bottom}.

It is crucial to note that the diffusive operator that contributes to ``\mathcal{L}^n_{ij}``
does not include fluxes across boundary faces.
In particular, ``\mathcal{L}^n_{ij}`` in \eqref{implicitoperatormatrix} enforces a
no-flux condition across the top and bottom faces.
Accordingly, fluxes through boundary faces due either to Dirichlet (Value)
boundary conditions or non-zero fluxes are accounted for
by adding the contribution of the flux diverence across
the top and bottom face to ``R \left ( \Phi \right )``.
