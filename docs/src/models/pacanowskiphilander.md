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

# The Pacanowski-Philander (1981) model

[Pacanowski and Philander (1981)][PP81] propose a simple one-dimensional model
for equatorial boundary layers dominated by mechanical turbulent mixing.
In their model, horizontal velocity, temperature, and salinity are governed by

```math
\beqs
U_t =   f V - \d_z \overline{w u} - F^u \c \\
V_t = - f U - \d_z \overline{w v} - F^v \c \\
S_t =       - \d_z \overline{w s} - F^S \c \\
T_t =       - \d_z \overline{w \theta} - F^T \c 
\eeqs
```

where uppercase variables are resolved, mean quantities, 
and lowercase variables are unresolved perturbations. A key variable in the
Pacaonwski and Philander formulation is the local Richardson number, defined by

```math
Ri = - \frac{g \rho_z}{\rho_0 \left ( U_z^2 + V_z^2 \right )} \c
```

where ``g`` is the gravitational constant, ``\rho_0`` is a reference density,
and ``\rho`` is density. For simplicity, we use a linear equation of state
between temperature, salinity, and density,

```math
\rho(z, t) = \rho_0 \left [ 
  \alpha \left ( T-T_0 \right ) + \beta \left ( S-S_0 \right ) \right ] \c
```
where typical values for the temperature and salinity expansion coefficients 
are ``\alpha = 2 \times 10^{-4} \r{K m^3 / kg}`` and ``\beta = 1``.

## Eddy diffusivities for momentum, temperature, and salinity

The core of the [PP81][PP81] model is to parameterize vertical turbulent fluxes 
with an eddy diffusivity/eddy viscosity, such that 

```math
\overline{w \phi} = \kappa_{\Phi} \d_z \Phi \c
```

where ``\kappa_{\Phi}`` is the eddy diffusivity or viscosity of the quantity 
``\Phi``.

The eddy viscosity for both ``U`` and ``V`` is

```math
\kappa_U = \nu_0 + \frac{\nu_1}{\left ( 1 + c Ri \right )^n} \c
```

while the eddy diffusivity for ``T`` and ``S`` are

```math
\kappa_T = \kappa_0 + \frac{ \kappa_1 }{ \left ( 1 + c Ri \right )^{n+1}} \p
```

This parameterization implies an ``Ri``-dependent turbulent Prandtl number of

```math
Pr = \frac{\kappa_U}{\kappa_T} \approx \left ( 1 + c Ri \right ) \c
```

when ``Ri`` is small. Typical values for the parameters (see [CV12][CV12]) are

* ``\nu_0 = 10^{-4} \r{cm/s}``,
* ``\nu_1 = 10^{-2} \r{cm/s}``,
* ``\kappa_0 = 10^{-5} \r{cm/s}``,
* ``\kappa_1 = 10^{-2} \r{cm/s}``,
* ``c = 5``, and
* ``n=2``.

[Pacanowski and Philander (1981)]: https://journals.ametsoc.org/doi/abs/10.1175/1520-0485(1981)011%3C1443:POVMIN%3E2.0.CO;2)
[PP81]: https://journals.ametsoc.org/doi/abs/10.1175/1520-0485(1981)011%3C1443:POVMIN%3E2.0.CO;2)
[CV12]: https://books.google.com/books?id=AAfoCAAAQBAJ
