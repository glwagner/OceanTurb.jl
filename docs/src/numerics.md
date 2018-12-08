```math
\newcommand{\c}{\, ,}

\newcommand{\r}[1]{\mathrm{#1}}

\newcommand{\ee}{\mathrm{e}}

\newcommand{\beq}{\begin{equation}}
\newcommand{\eeq}{\end{equation}}

\newcommand{\beqs}{\begin{gather}}
\newcommand{\eeqs}{\end{gather}}
```


# Numerical modeling of the oceanic boundary layer

## Basic form

Models for the oceanic boundary layer are partial differential equations of 
the form

```math
\beq
\phi_t = L \phi + N(\phi) \c
\label{mathematicalform}
\eeq
```

where ``phi`` is a variable like velocity, temperature, or salinity, ``L`` is 
a linear operator, and ``N`` is a nonlinear operator.

## Time-stepping

An explicit forward Euler time integration scheme discretizes
\eqref{mathematicalform} in time with

```math
\beq
\phi^{n+1} = \phi^{n} + \Delta t \left [ L \phi^n + N(\phi^n) \right ] \,
\eeq 
```

where the superscripts ``n`` and ``n+1`` denote the solution at 
time-step ``n`` and ``n+1``, respectively.
