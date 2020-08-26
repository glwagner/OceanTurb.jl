#
# Generic diffusivity functions
#

@inline KU(m, i) = KU₀(m) + Cᴷu(m, i) * onface(m.state.K, i)
@inline KV(m, i) = KV₀(m) + Cᴷv(m, i) * onface(m.state.K, i)
@inline KT(m, i) = KT₀(m) + CᴷT(m, i) * onface(m.state.K, i)
@inline KS(m, i) = KS₀(m) + CᴷS(m, i) * onface(m.state.K, i)
@inline Ke(m, i) = Ke₀(m) + Cᴷe(m, i) * onface(m.state.K, i)

#
# Background diffusivity parameters
#

Base.@kwdef struct BackgroundDiffusivities{T} <: AbstractParameters
    KU₀ :: T = 1e-6 # Background viscosity for momentum [m s⁻²]
    KC₀ :: T = 1e-6 # Background diffusivity for tracers and TKE [m s⁻²]
end

@inline KU₀(m) = m.background_diffusivities.KU₀
@inline KV₀(m) = m.background_diffusivities.KU₀
@inline KC₀(m) = m.background_diffusivities.KC₀
@inline KT₀(m) = m.background_diffusivities.KC₀
@inline KS₀(m) = m.background_diffusivities.KC₀
@inline Ke₀(m) = m.background_diffusivities.KC₀

#
# Eddy diffusivity parameters
#

"""
    struct SinglePrandtlDiffusivities{T} <: AbstractParameters

A `simple' Prandtl number model that assumes the TKE diffusivity is identical
to the tracer diffusivity. The momentum diffusivity parameter is Cᴷu, and the 
tracer diffusivities equal to the momentum diffusivity divided by CᴷPr.
"""
Base.@kwdef struct SinglePrandtlDiffusivities{T} <: AbstractParameters
     Cᴷu :: T = 0.1   # Diffusivity parameter for velocity
    CᴷPr :: T = 0.74  # Constant Prandtl number for tracers and TKE
end

const SPD = SinglePrandtlDiffusivities

@inline Cᴷu(m::Model{L, <:SPD}, i) where L = m.eddy_diffusivities.Cᴷu
@inline Cᴷv(m::Model{L, <:SPD}, i) where L = m.eddy_diffusivities.Cᴷu
@inline Cᴷc(m::Model{L, <:SPD}, i) where L = m.eddy_diffusivities.Cᴷu / m.eddy_diffusivities.CᴷPr
@inline CᴷT(m::Model{L, <:SPD}, i) where L = Cᴷc(m, i)
@inline CᴷS(m::Model{L, <:SPD}, i) where L = Cᴷc(m, i)
@inline Cᴷe(m::Model{L, <:SPD}, i) where L = Cᴷc(m, i)


"""
    struct IndependentDiffusivities{T} <: AbstractParameters

A diffusivity model in which momentum, tracers, and TKE
each have their own diffusivity parameter.
"""
Base.@kwdef struct IndependentDiffusivities{T} <: AbstractParameters
     Cᴷu :: T = 0.0274 # Diffusivity parameter for velocity
     Cᴷc :: T = 0.0498 # Diffusivity parameter for tracers
     Cᴷe :: T = 0.0329 # Diffusivity parameter for TKE
end

const IDP = IndependentDiffusivities 

@inline Cᴷu(m::Model{L, <:IDP}, i) where L = m.eddy_diffusivities.Cᴷu
@inline Cᴷv(m::Model{L, <:IDP}, i) where L = m.eddy_diffusivities.Cᴷu
@inline Cᴷc(m::Model{L, <:IDP}, i) where L = m.eddy_diffusivities.Cᴷc
@inline CᴷT(m::Model{L, <:IDP}, i) where L = Cᴷc(m, i)
@inline CᴷS(m::Model{L, <:IDP}, i) where L = Cᴷc(m, i)
@inline Cᴷe(m::Model{L, <:IDP}, i) where L = m.eddy_diffusivities.Cᴷe


"""
    struct RiDependentDiffusivities{T} <: AbstractParameters

A diffusivity model in which momentum, tracers, and TKE
each have independent Richardson number dependent diffusivities.

The Richardson number is

    ``Ri = ∂z B / ( (∂z U)² + (∂z V)² )`` ,

where ``B`` is buoyancy and ``∂z`` denotes a vertical derviative.

The Richardson-number dependent diffusivities are multiplied by the stability
function

    ``σ(Ri) = σ⁰ + σᵟ * step(Ri, Riᶜ, Riʷ)``
    
where ``σ⁰``, ``σᵟ``, ``Riᶜ``, and ``Riʷ`` are free parameters,
and ``step`` is a smooth step function defined by

    ``step(x, c, w) = 1/2 * (1 + tanh((x - c) / w))``.
"""
Base.@kwdef struct RiDependentDiffusivities{T} <: AbstractParameters
     Cᴷu⁻  :: T = 0.02   # Convection diffusivity parameter for velocity
     Cᴷuᵟ  :: T = 0.01   # Shift diffusivity parameter for velocity
     Cᴷc⁻  :: T = 0.04   # Convection diffusivity parameter for tracers
     Cᴷcᵟ  :: T = 0.01   # Shift diffusivity parameter for tracers
     Cᴷe⁻  :: T = 0.02   # Convection diffusivity parameter for TKE
     Cᴷeᵟ  :: T = 0.01   # Shift diffusivity parameter for TKE
     CᴷRiʷ :: T = 0.1    # Ri width parameter
     CᴷRiᶜ :: T = 0.1    # "Central" Ri parameter
end

const RiD = RiDependentDiffusivities

@inline function Richardson_number(m, i)
    N² = ∂B∂z(m, i)
    return ifelse(N² == 0, 0.0, N² / shear_squared(m, i))
end

@inline step(x, c, w) = 1/2 * (1 + tanh((x - c) / w))

@inline stability_function(Ri, σ⁻, σᵟ, c, w) = σ⁻ + σᵟ * step(Ri, c, w)

@inline Cᴷu(m::Model{L, <:RiD}, i) where L = stability_function(
                                                                Richardson_number(m, i),
                                                                m.eddy_diffusivities.Cᴷu⁻,
                                                                m.eddy_diffusivities.Cᴷuᵟ,
                                                                m.eddy_diffusivities.CᴷRiᶜ,
                                                                m.eddy_diffusivities.CᴷRiʷ,
                                                               )

@inline Cᴷv(m::Model{L, <:RiD}, i) where L = Cᴷu(m, i)

@inline Cᴷc(m::Model{L, <:RiD}, i) where L = stability_function(
                                                                Richardson_number(m, i),
                                                                m.eddy_diffusivities.Cᴷc⁻,
                                                                m.eddy_diffusivities.Cᴷcᵟ,
                                                                m.eddy_diffusivities.CᴷRiᶜ,
                                                                m.eddy_diffusivities.CᴷRiʷ,
                                                               )

@inline CᴷT(m::Model{L, <:RiD}, i) where L = Cᴷc(m, i)
@inline CᴷS(m::Model{L, <:RiD}, i) where L = Cᴷc(m, i)

@inline Cᴷe(m::Model{L, <:RiD}, i) where L = stability_function(
                                                                Richardson_number(m, i),
                                                                m.eddy_diffusivities.Cᴷe⁻,
                                                                m.eddy_diffusivities.Cᴷeᵟ,
                                                                m.eddy_diffusivities.CᴷRiᶜ,
                                                                m.eddy_diffusivities.CᴷRiʷ,
                                                               )
