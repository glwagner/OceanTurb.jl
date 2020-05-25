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


