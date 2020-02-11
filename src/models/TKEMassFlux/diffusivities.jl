#
# Background diffusivity parameters
#

Base.@kwdef struct BackgroundDiffusivities{T} <: AbstractParameters
    KU₀ :: T = 1e-6 # Background viscosity for momentum
    KC₀ :: T = 1e-6 # Background diffusivity for tracers
end

@inline KU₀(m) = m.background_diffusivities.KU₀
@inline KC₀(m) = m.background_diffusivities.KC₀

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
    CᴷPr :: T = 0.74  # Diffusivity parameter for velocity
end

const SPD = SinglePrandtlDiffusivities

@inline KU(m::Model{L, <:SPD}, i) where L = KU₀(m) + m.eddy_diffusivities.Cᴷu * onface(m.state.K, i)
@inline KV(m::Model{L, <:SPD}, i) where L = KU(m, i)

# All tracers share diffusivity.
@inline KC(m::Model{L, <:SPD}, i) where L = KC₀(m) + m.eddy_diffusivities.Cᴷu /
                                                     m.eddy_diffusivities.CᴷPr * onface(m.state.K, i)

@inline KT(m::Model{L, <:SPD}, i) where L = KC(m, i)
@inline KS(m::Model{L, <:SPD}, i) where L = KC(m, i)
@inline Ke(m::Model{L, <:SPD}, i) where L = KC(m, i)


"""
    struct IndependentDiffusivities{T} <: AbstractParameters

A diffusivity model in which momentum, tracers, and TKE
each have their own diffusivity parameter.
"""
Base.@kwdef struct IndependentDiffusivities{T} <: AbstractParameters
     Cᴷu :: T = 0.1   # Diffusivity parameter for velocity
     Cᴷc :: T = 0.15  # Diffusivity parameter for tracers
     Cᴷe :: T = 0.15  # Diffusivity parameter for TKE
end

const IDP = IndependentDiffusivities 

# Momentum diffusivities:
@inline KU(m::Model{L, <:IDP}, i) where L = KU₀(m) + m.eddy_diffusivities.Cᴷu * onface(m.state.K, i)
@inline KV(m::Model{L, <:IDP}, i) where L = KU(m, i)

# TKE diffusivity
@inline Ke(m::Model{L, <:IDP}, i) where L = KC₀(m) + m.eddy_diffusivities.Cᴷe * onface(m.state.K, i)

# Tracer diffusivities:
@inline KC(m::Model{L, <:IDP}, i) where L = KC₀(m) + m.eddy_diffusivities.Cᴷc * onface(m.state.K, i)

@inline KT(m::Model{L, <:SPD}, i) where L = KC(m, i)
@inline KS(m::Model{L, <:SPD}, i) where L = KC(m, i)
