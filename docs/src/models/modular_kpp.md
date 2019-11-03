# Modular K-Profile Parameterization

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
\newcommand{\Ek}        {\mathrm{Ek}}
\newcommand{\SL}        {\mathrm{SL}}
\newcommand{\K}         {\mathcal{E}}
\newcommand{\W}         {\mathcal{W}}

\newcommand{\btau}      {\b{\tau}} % wind stress vector

% Model functions and constants
\renewcommand{\F}[2]    {\Upsilon^{#1}_{#2}}
\renewcommand{\C}[2]    {C^{#1}_{#2}}

\newcommand{\uwind}     {u_\star}
\newcommand{\ubuoy}     {w_\star}

\newcommand{\NL}        {NL}
```

In the `ModularKPP` module, horizontally-averaged vertical turbulent
fluxes are modeled with the combination of a local diffusive flux and a non-local
non-diffusive flux:


```math
\beq
\overline{w \phi} = - K_\Phi \d_z \Phi + \NL_\Phi \c
\eeq
```

where the depth dependence of the eddy diffusivity ``K_\Phi`` is 

```math
\beq
K_\Phi \propto h \, \W_\Phi \, \F{d}{}(d) \, ,
\eeq
```

where ``\W_\Phi`` is a turbulent velocity scale that in general depends on 
``\Phi``, the quantity being diffused, ``d \equiv - z / h`` ,
and ``h`` is the 'mixing layer depth'.
Typically ``\F{d}{}`` is the cubic polynomial

```math
\beq
\F{d}{}(d) = d ( 1 - d )^2 \, ,
\eeq
```

however, `ModularKPP` permits experimentation with different forms.


The formulation of diffusivity as the product of a magnitude with a with a shape or 
'profile' function gives rise to the name.
``K``-profile parameterization.

The non-local flux term ``\NL_\Phi`` models the effects of convective
plumes.

``K``-profile schemes with a non-local flux term thus have three basic components:

1. A model for the mixing layer depth ``h``, over which ``K_\Phi > 0``.
2. A model for the magnitude of the diffusivity, ``K``
3. A model or "shape function" that determines the dependence of ``K`` as a function of ``d=-z/h``.
4. A model for the non-local flux, ``\NL_\Phi``.

# Model instantiaton

A `ModularKPP.Model` is instantiated in the default configuration by
writing 

```julia
using OceanTurb
model = ModularKPP.Model()
```

or,

```julia
using OceanTurb

model = Model(grid = UniformGrid(N=10, L=1.0),
         constants = Constants(),
       diffusivity = LMDDiffusivity(),
      nonlocalflux = LMDCounterGradientFlux(),
       mixingdepth = LMDMixingDepth(),
          kprofile = StandardCubicPolynomial(),
           stepper = :BackwardEuler,
               bcs = ModelBoundaryConditions(eltype(grid)),
           forcing = Forcing())
```

This builds a model on a `UniformGrid` with `N=10` grid points and `L=1.0` meters deep.
The keyword arguments `diffusivity`, `nonlocalflux`, `mixingdepth`, and `kprofile` correspond 
to a specific ``K``-profile configuration proposed by
[Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872):

* `diffusivity = LMDDiffusivity()` determines the turbulent velocity scale ``W_\Phi`` 
    using the prescription proposed by LMD94

* `nonlocalflux = LMDCounterGradientFlux()` determines the nonlocal flux ``NL_\Phi``
    using the prescription proposed by LMD94

* `mixingdepth = LMDMixingDepth()` determines the mixing depth ``h`` using the bulk
    Richardson number criterion proposed by LMD94

* `kprofile = StandardCubicPolynomial()` uses the cubic polynomial ``d(1-d)^2`` to set the primary 
    depth dependence of ``K_\Phi``, as proposed by LMD94.

More subcomponent choices and details about their consequences are described in 
[Sub-components of `ModularKPP.Model`](@ref).

The keyword arguments `constants`, `stepper`, `bcs`, and `forcing` configure the model constants
time stepper, boundary conditions, and forcing function.
The only useful time-stepper at the moment is `:BackwardEuler`.
The procedures for setting boundary conditions and defining forcing functions are described in 
[Setting boundary conditions](@ref) and [Defining forcing functions](@ref).

The default set of constants is returned by `constants = Constants()`, or

```julia
constants = Constants(
     α = 2.5e-4, # thermal expansion coefficient [C⁻¹]
     β = 8e-5,   # haline contraction coefficient [psu⁻¹]
    ρ₀ = 1035,   # reference density [kg m⁻³]
    cP = 3992,   # heat capacity `cP` [...]
     f = 0,      # Coriolis parameter `f` [s⁻¹]
     g = 9.81    # gravitational acceleration `g` [m² s⁻¹]
)
```

# Setting boundary conditions

Two basic methods may be used to set boundary conditions. 

## Constant and standard boundary conditions

For boundary conditions consisting of constant surface fluxes or constant
bottom gradients, 

```julia
using OceanTurb

model = ModularKPP.Model()

model.bcs.U.top = BoundaryCondition(Flux, -1e-4)
model.bcs.T.top = BoundaryCondition(Flux, 1e-4)
model.bcs.T.bottom = BoundaryCondition(Gradient, model.constants.α * model.constants.g * 1e-5)
```

will, for example, set the top boundary condition on temperature `T` to a positive flux
of ``10^{-4} \, \mathrm{m K \, s^{-2}}``, a bottom temperature gradient that corresponds
to a bottom buoyancy gradient of ``N^2 = 10^{-5} \, \mathrm{s^{-2}}``,
and a top boundary condition on the horizontal velocity `U` to a negative flux.

In an oceanic scenario, a positive surface temperature flux of ``Q_\theta = 10^{-4}`` 
is strongly destabilizing,
corresponding to a heat flux of ``Q_h = \rho_0 c_P Q_\theta \approx 413 \, \mathrm{W \, m^{-2}}``,
or in ordinary oceanographic parlance a 'heating' of ``-413 \, \mathrm{W \, m^{-2}}``.
(A positive surface flux extracts a quantity from the oceanic domain below; therefore 
positive temperature flux acts to cool and destabilize at the ocean surface.
This convention is standard --- an upward velocity leads to a positive flux, for example ---
but is opposite the ordinary convention in oceanography.)

A negative flux of velocity accelerates surface fluid in the positive ``x``-direction.
A velocity flux, or kinematic stress of ``Q_u = -10^{-4}`` corresponds to a friction velocity
of ``u_\star = | \boldsymbol{Q}_u |^{1/2} = 0.01 \, \mathrm{m \, s^{-1}}`` and a dynamic stress of
``\boldsymbol{\tau} = \rho_0 \boldsymbol{Q}_u \approx -10^{-1} \, \mathrm{N \, m^{-2}}``.

## More complex boundary conditions

For non-standard or more complicated boundary conditions that are enforced, for example, by
time-dependent or nonlinear functions, a variable's boundary condition must be generated
prior to model instantiation and passed to the model constructor.
To set a time-dependent surface flux of temperature for example, write

```julia
using OceanTurb

# Functions-as-boundary-conditions take a single argument of type `ModularKPP.Model`.
fun_flux(model) = 1e-8 * cos(2π/day * model.clock.time)

# Wrap `fun_flux` in a `BoundaryCondition` and specify its application as a flux.
top_temperature_bc = BoundaryCondition(Flux, fun_flux)

# Instantiate boundary conditions for temperature with the flux function on top.
temperature_bcs = BoundaryConditions(top=top_temperature_bc)

# Instantiate a model with the indicated temperature boundary condition and default
# boundary conditions for all other variables.
model = Model(bcs = ModelBoundaryConditions(T=temperature_bcs))

# Constant boundary conditions of default type on other variables are still settable.
model.bcs.T.bottom = BoundaryCondition(Gradient, model.constants.α * model.constants.g * 1e-5)
```

# Defining forcing functions

Forcing functions have the signature `forcing_func(model, i)`, where `model::ModularKPP.Model`, and `i`
is the grid point at which the forcing is applied.
For example, to apply a body force on `U`, write

```julia
using OceanTurb

@inline body_force(model, i) = @inbounds -1e-1 * model.grid.zc[i] / model.grid.L

model = Model(forcing = Forcing(U=body_force))
```

This instantiates a model with the specified body force applied to ``U``, such that the ``U`` equation becomes

```math
\beq
\partial_t U - f V = \partial_z \left ( K_U \partial_z U \right ) - 10^{-1} \frac{z}{L} \, .
\eeq
```

The annotations `@inline` tells the julia compiler to "inline" the function, which typically 
increases performance, and the `@inbounds` annotation instructs the compiler to elide 
bounds checking when indexing the range `model.grid.zc`, which also saves time provided
that `body_force` is never called with `i` out of bounds.

# Sub-components of `ModularKPP.Model`

## Mixing depth models

The mixing depth model is configured via the keyword argument `mixingdepth` in the `ModularKPP.Model`
constructor.

### CVMix mixing depth model

The CVMix mixing depth model is instiated by writing

```julia
mixingdepth = LMDMixingDepth()
```

The
[CVMix](https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf)
mixing depth model uses the 'bulk Richardson number' criterion proposed by
[Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872).
This model is described in [Mixing depth model in CVMix KPP](@ref).

### ROMS mixing depth model

The ROMS mixing depth model is instantiated by writing

```julia
mixingdepth = ROMSMixingDepth()
```

The mixing depth model used by the [Regional Ocean Modeling System (ROMS)](https://www.myroms.org)
is described in appendix B of
[McWilliams et al (2009)](https://journals.ametsoc.org/doi/full/10.1175/2009JPO4130.1).
The model introduces a 'mixing function' ``\mathbb{M}``, which is increased
by shear and convection and decreased by stable stratification and rotation.
``\mathbb{M}`` is defined

```math
\beq \label{stabilization}
\mathbb M(z) = \int_z^0 \F{\SL}{}(z') \left [
      \left ( \d_z \b{U} \right )^2 - \frac{\d_z B}{\C{\Ri}{}} - \C{\Ek}{} f^2
    \right ] \, \mathrm{d} z'
    - \C{\K}{} \ubuoy^\dagger N^\dagger \c
\eeq
```
where

```math
\beq
\ubuoy^\dagger(z) \equiv \max \left (0, -z F_b \right )^{1/3} \c
  \quad \mathrm{and} \quad
N^\dagger(z) \equiv \max \left (0, \d_z B \right )^{1/2} \p
\eeq
```
Typically, the mixing function ``\mathbb M(z)`` increases from 0 at ``z=0``
into the well-mixed region immediately below the surface due to
``\left ( \d_z \b{U} \right )^2`` and ``\ubuoy^\dagger N^\dagger``
during convection, and decreases to negative values in the stratified
region below the mixing layer due to the stabilizing action of ``-\d_z B / \C{\Ri}{}``.
The boundary layer depth is defined as the first nonzero depth where ``\mathbb M(z) = 0``.

Finally, the 'surface layer exclusion' function,

```math
\beq \label{exclusion}
\F{\SL}{} \equiv - \frac{z}{\C{\SL}{} h - z} \c
\eeq
```
acts to exclude the values of
``\left ( \d_z \b{U} \right )^2`` and ``\d_z B`` at the top of the boundary
layer from influencing the diagnosed boundary layer depth.

[McWilliams et al (2009)](https://journals.ametsoc.org/doi/full/10.1175/2009JPO4130.1)
suggest
``\C{\SL}{} = 0.1``, ``\C{\K}{} = 5.07``, ``\C{\Ri}{} = 0.3``, and ``\C{\Ek}{} = 211``
for the free parameters in \eqref{stabilization}.

## Diffusivity models

### [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)

The diffusivity model proposed by 
[Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
(LMD94) is instantiated by writing

```julia
diffusivity = LMDDiffusivity()
```

The LMD94 diffusivity model prescribes the turbulent velocity scale
``\W_\Phi(d)`` in the generic ``K``-profile formulation,

```math
\beq
K_\Phi \propto h \W^\text{LMD94}_\Phi(d) \F{d}{}(d) \, .
\eeq
```

The formulation of ``\W^\text{LMD94}_\Phi(d)`` is described in 
[``K``-Profile model in CVMix KPP](@ref).

`ModularKPP.Model` permits a range of shape functions ``\F{d}{}(d)`` to be used with the
LMD94 turbulent velocity scale ``\W_\Phi(d)``.

### Holtslag (1998)

The diffusivity model proposed by Holtslag in 1998 and described in 
[Siebesma et al (2007)](https://journals.ametsoc.org/doi/full/10.1175/JAS3888.1)
is instantiated by writing

```julia
diffusivity = HoltslagDiffusivity()
```

The Holtslag diffusivity uses the simple turublent velocity scale,

```math
\beq
\W^\text{Holtslag} = \C{\tau}{} \ubuoy \left [ \left ( \frac{\uwind}{\ubuoy} \right )^3 
    + \C{\tau b}{} d \right ]^{1/3} \, .
\eeq
```

[Siebesma et al (2007)](https://journals.ametsoc.org/doi/full/10.1175/JAS3888.1), which
pair the turbulent velocity scale ``\W^\text{Holtslag}`` with a cubic shape function, 
a diagnostic plume model and simple mixing depth model, suggest
``\C{\tau}{} = 0.4`` and ``\C{\tau b}{} = 15.6``.

## Non-local flux models

### 'Countergradient flux' model

The counter gradient flux model proposed by 
[Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
is instantiated by writing

```julia
nonlocalflux = LMDCounterGradientFlux()
```

As described in 
['Countergradient' non-local flux model in CVMix KPP](@ref),
the non-local countergradient flux is defined only for ``T`` and ``S``, and is

```math
\beq
NL_\phi = \C{\NL}{} F_\phi d (1 - d)^2 \c
\eeq
```
where ``d = -z/h`` is a non-dimensional depth coordinate and ``\C{\NL}{} = 6.33``.

### Diagnostic plume model

The diagnostic plume model proposed by
[Siebesma et al (2007)](https://journals.ametsoc.org/doi/full/10.1175/JAS3888.1)
is instantiated by writing

```julia
nonlocalflux = DiagnosticPlumeModel()
```

The diagnostic plume model
integrates equations that describe the quasi-equilibrium vertical momentum and 
tracer budgets for plumes that plunge downwards from the ocean surface
due to destabilizing buoyancy flux. 

In the diagnostic plume model, the non-local flux of a tracer ``\Phi`` is parameterized as
```math
\beq
    NL_\phi = \C{a}{} \breve W \left ( \Phi - \breve \Phi \right ) \, ,
\eeq
```
where ``\C{a}{} = 0.1`` is the plume area fraction, ``\breve W`` is the plume
vertical velocity, and ``\breve \Phi`` is plume-averaged concentration of the 
tracer ``\phi``.
When using a plume model in `ModularKPP`, ``\Phi`` must be interpreted as
the average concentration of ``\phi`` in the environment, excluding plume regions.
in `OceanTurb.jl`'s implementation of the
[Siebesma et al (2007)](https://journals.ametsoc.org/doi/full/10.1175/JAS3888.1)
the plume vertical velocity variance ``\breve W^2`` is used as a diagnostic variable, 
rather than ``\breve W``. Due to this, the nonlocal flux in `OceanTurb.jl` becomes
```math
\beq
    NL_\phi = -\C{a}{} \sqrt{\breve W^2} \left ( \Phi - \breve \Phi \right ) \, ,
\eeq
```
where we assume that ``\breve W \le 0``.


#### Continuous plume equations

The diagnostic, steady-state plume-averaged temperature and salinity budgets boil down to
```math
\begin{gather}
    \d_z \breve T = - \epsilon \left ( \breve T - T \right ) \, , \\
    \d_z \breve S = - \epsilon \left ( \breve S - S \right ) \, ,
\end{gather}
```
where ``\breve T`` and ``\breve S`` are the plume-averaged temperature and salinity, 
and ``T`` and ``S`` are the environment-averaged temperature and salinity and 
``\epsilon(z, h)`` is the parameterized entrainment rate, 
```math
\beq
\epsilon(z) = \C{\epsilon}{} 
    \left [ \frac{1}{\Delta c_N - z} + \frac{1}{\Delta c_N + \left ( z + h \right )} \right ] \, ,
\eeq
```
where ``\C{\epsilon}{} = 0.4``, ``\Delta c_N`` is the spacing between the boundary and the topmost
cell interface, and ``h`` is the mixing depth determined via the chosen mixing depth model.

The budget for plume vertical momentum is

```math
\beq
    \d_z \breve W^2 = \C{b}{w} \left ( \breve B - B \right ) - \C{\epsilon}{w} \epsilon \, \breve W^2
\eeq
```
where ``\breve B = \alpha \breve T - \beta \breve S`` is the plume-averaged buoyancy and 
``B = \alpha T - \beta S`` is the environment-averaged buoyancy, ``\C{b}{w} = 2.86``, 
and ``\C{\epsilon}{w} = 0.572``.

#### Surface layer plume initialization model and numerical implementation

The plume equations require boundary conditions at ``z=0``, or an initialization model 
at the topmost node in the interior of the domain. 
Note that ``\breve T`` and ``\breve S`` are defined at cell centers, and 
``\breve W^2`` is defined at cell interfaces.

The plume-averaged tracer concentration in the topmost cell is parameterized in terms
of the tracer flux across the top boundary, ``Q_\phi``, with the formula
```math
\beq
    \breve \Phi(z=z_N) = \Phi(z=z_N) - \C{\alpha}{} \frac{Q_\phi}{\sigma_w(z_N)} \, ,
\eeq
```
where ``\C{\alpha}{} = 1.0``, and ``\sigma_w(z)`` is an empirical expression for the 
vertical velocity standard deviation,
```math
\beq
    \sigma_w = \left ( \C{\sigma \tau}{} \uwind^3 + \C{\sigma b}{} \ubuoy^3 d \right )^{1/3} \left ( 1 - d \right )^{1/2} \, ,
\eeq
```
with ``\C{\sigma \tau}{} = 2.2`` and ``\C{\sigma b}{} = 1.32``.

The boundary condition on plume vertical momentum prescribes no penetration
through the ocean surface,
```math
\beq
    \breve W(z=0) = 0 \, .
\eeq
```

The advetion terms in the diagnostic plume equations are discretized with an downwind 
scheme, which permits integration of each equation from the surface downward.
The plume temperature advection term, for example, is defined at cell centers
and discretized with

```math
\beq
\left ( \partial_z \breve T \right )_{i+1} = \frac{\breve T_{i+1} - \breve T_i}{\Delta c_{i+1}} \, .
\eeq
```

At the same time, the terms on the right side of the plume temperature conservation 
equation are evaluated at cell center `i+1`, which leads to the integral

```math
\beq
    \breve T_i = \breve T_{i+1} + \Delta c_{i+1} \epsilon \left ( z_{c, i+1} \right ) 
        \left ( \breve T_{i+1} - T_{i+1} \right ) \, ,
\eeq
```

that determines the plume temperature at cell center `i` with respect to plume temperature,
environment temperature, and entrainment quantities evaluated at cell center `i+1`.
The plume salinity conservaiton equation is discretized analogously.

The downwind discretization of the plume vertical momentum advection term is

```math
\beq
\left ( \partial_z \breve W^2 \right )_{i+1} = \frac{ \breve W^2_{i+1} - \breve W^2_i }{\Delta f_{i+1}} \, .
\eeq
```

The plume vertical momentum is defined at cell interfaces. This means that discretizing the
right side of the plume vertical momentum equation requires interpolating the buoyancy
field from cell centers to cell interfaces.
The plume vertical velocity equation is thus

```math
\beq
    \breve W^2_i = \breve W^2_{i+1} - \Delta f_{i+1} \C{b}{w} \frac{1}{2} \left ( \breve B_i + \breve B_{i+1} - B_i - B_{i+1} \right ) + \Delta f_{i+1} \C{\epsilon}{w} \epsilon \left ( z_{f, i+1} \right ) \breve W^2_{i+1} \, .
\eeq
```

The plume integration is stopped at grid point `i` when ``\breve W^2_{i+1} < 0``.

To numerically integrate the environment-averaged tracer conservation equation, 
the mass flux term is divided into two components,
```math
\beq
\partial_t \Phi - \partial_z \left ( K_\Phi \partial_z \Phi \right ) 
                + \partial_z \left ( -\C{a}{} \sqrt{\breve W^2} \Phi \right )
                = - \C{a}{} \sqrt{\breve W^2} \breve \Phi \, ,
\eeq
```
where the diffusivity and mass flux term on the left are integrated implicitly in time, and the
mass flux term on the right is integrated explicitly in time.
