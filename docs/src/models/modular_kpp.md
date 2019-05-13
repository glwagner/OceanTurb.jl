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

\newcommand{\btau}      {\b{\tau}} % wind stress vector

% Model functions and constants
\renewcommand{\F}[2]      {\Upsilon^{#1}_{#2}}
\renewcommand{\C}[2]      {C^{#1}_{#2}}

\newcommand{\uwind}     {\omega_{\tau}}
\newcommand{\ubuoy}     {\omega_b}

\newcommand{\NL}        {NL}
```

In `ModularKPP` module, horizontally-averaged vertical turbulent
fluxes are modeled with the combination of a local diffusive flux and a non-local
non-diffusive flux:


```math
\beq
\overline{w \phi} = - K_\Phi \d_z \Phi + \NL_\Phi \c
\eeq
```

where the depth dependence of the eddy diffusivity ``K_\Phi`` is modeled
with a shape or 'profile' function, which gives rise to the name
``K``-profile parameterization.
The non-local flux term ``\NL_\Phi`` models the effects of convective
plumes.

``K``-profile schemes with a non-local flux term thus have three basic components:

1. A model for the mixing layer depth ``h``, over which ``K_\Phi > 0``.
2. A model for the diffusivity ``K``, which includes the ``K``-profile as well as its magnitude.
3. A model for the non-local flux, ``\NL_\Phi``.

## Mixing depth models

### CVMix mixing depth model

The
[CVMix](https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf)
mixing depth model uses the 'bulk Richardson number' criterion proposed by
[Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872).
This model is described in [Mixing depth model in CVMix KPP][@ref].

### ROMS mixing depth model

The mixing depth model used by the [Regional Ocean Modeling System (ROMS)](https://www.myroms.org)
is described in appendix B of
[McWilliams et al (2009)](https://journals.ametsoc.org/doi/full/10.1175/2009JPO4130.1).
The model introduces a 'mixing function' ``\mathbb{M}``, which is increased
by shear and convection and decreased by stable stratification and rotation.
The stabilization function is defined as

```math
\beq \label{stabilization}
\mathbb M(z) = \int_z^0 \F{\SL}{}(z') \left [
      \left ( \d_z \b{U} \right )^2 - \frac{\d_z B}{\C{\Ri}{}} - \C{\Ek}{} f^2
    \right ] \, \mathrm{d} z'
    - \C{\K}{} \omega_b^\dagger N^\dagger \c
\eeq
```
where

```math
\beq
\omega_b^\dagger(z) \equiv \max \left (0, -z F_b \right )^{1/3} \c
  \quad \mathrm{and} \quad
N^\dagger(z) \equiv \max \left (0, \d_z B \right )^{1/2} \p
\eeq
```
Typically, the mixing function ``\mathbb M(z)`` increases from 0 at ``z=0``
into the well-mixed region immediately below the surface due to
``\left ( \d_z \b{U} \right )^2`` and ``\omega_b^\dagger N^\dagger``
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
suggests
``\C{\SL}{} = 0.1``, ``\C{\K}{} = 5.07``, ``\C{\Ri}{} = 0.3``, and ``\C{\Ek}{} = 211``
for the free parameters in \eqref{stabilization}.

## Diffusivity models

## Non-local flux models

* [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872)
* [CVMix documentation](https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf)
* [McWilliams et al (2009)](https://journals.ametsoc.org/doi/full/10.1175/2009JPO4130.1)
* [Regional Ocean Modeling Systeml (ROMS)][https://www.myroms.org]

[^LMD94]: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872
[^MHS09]: https://journals.ametsoc.org/doi/full/10.1175/2009JPO4130.1
[^CVMix]: https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf
[^ROMS]: https://www.myroms.org
