# Pacanowski-Philander
```math
\newcommand{\c}     {\, ,}
\newcommand{\p}     {\, .}
\newcommand{\d}     {\partial}
\newcommand{\r}[1]  {\mathrm{#1}}
\newcommand{\ee}    {\mathrm{e}}
\newcommand{\beq}   {\begin{equation}}
\newcommand{\eeq}   {\end{equation}}
\newcommand{\beqs}  {\begin{gather}}
\newcommand{\eeqs}  {\end{gather}}
\newcommand{\Ri}    {\mathrm{Ri}}
```
In the model proposed by
[Pacanowski and Philander (1981)](https://journals.ametsoc.org/doi/abs/10.1175/1520-0485(1981)011%3C1443:POVMIN%3E2.0.CO;2),
turbulent fluxes are diffusive, so that

```math
\overline{w \phi} = K_\Phi \d_z \Phi \c
```

where the diffusivity for velocity fields, ``K_U``, is

```math
\beq \label{momentumdiffusivity}
K_U = \nu_0 + \frac{\nu_1}{\left ( 1 + c \Ri \right )^n} \c
\eeq
```

while the diffusivity for tracer fields is

```math
\beq \label{tracerdiffusivity}
K_T = \kappa_0 + \frac{\kappa_1}{\left ( 1 + c \Ri \right )^{n+1}} \p
\eeq
```

In \eqref{momentumdiffusivity} and \eqref{tracerdiffusivity}, the local Richardson number
``\Ri`` is defined

```math
\beq
Ri = \frac{\d_z B}{\left ( \d_z U \right )^2 + \left ( \d_z V \right )^2} \c
\eeq
```

in terms of the buoyancy ``B = - g \rho' / \rho_0``, where ``g`` is gravitational acceleration, ``\rho_0`` is a reference density, and ``\rho'`` is the density deviation therefrom.
With the linear equation of state

```math
\beq
\rho = \rho_0 \left [ 1 - \alpha \left ( T - T_0 \right ) + \beta \left ( S - S_0 \right ) \right ]
\eeq
```

near some reference temperature ``T_0`` and reference salinity ``S_0``, buoyancy ``B``
is given by

```math
\beq
B = g \left [ \alpha \left ( T - T_0 \right ) - \beta \left ( S - S_0 \right ) \right ] \c
\eeq
```

and its vertical derivative is

```math
\beq
\d_z B = g \left ( \alpha \d_z T - \beta \d_z S \right ) \p
\eeq
```

## Parameters

Typical values for the model parameters in PP
(see the text following equation 19 in chapter 3 of
[CV12](https://books.google.com/books?id=AAfoCAAAQBAJ))
are

|   Parameter   | Value         | Units                  |
|   :-------:   | :---:         | -----                  |
| ``\nu_0``     | ``10^{-4}``   | ``\r{m^2 \, s^{-1}}`` |
| ``\nu_1``     | ``10^{-2}``   | ``\r{m^2 \, s^{-1}}`` |
| ``\kappa_0``  | ``10^{-5}``   | ``\r{m^2 \, s^{-1}}`` |
| ``\kappa_1``  | ``10^{-2}``   | ``\r{m^2 \, s^{-1}}`` |
| ``c``         | ``5``         | none |
| ``n``         | ``2``         | none |
