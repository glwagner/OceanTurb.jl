# Fallbacks
@inline diffusivity_mixing_length(m, i) = mixing_length(m, i)
@inline dissipation_length(m, i) = mixing_length(m, i)

@inline shear_squared(m, i) = ∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2

#
# A "simple" mixing length
#

Base.@kwdef struct SimpleMixingLength{T} <: AbstractParameters
    Cᴸᵏ :: T = 0.4
    Cᴸᵇ :: T = 0.1
    Cᴸᵟ :: T = 1.0
end

@inline function mixing_length(m::Model{<:SimpleMixingLength}, i)
    # Two mixing lengths based on stratification and distance from surface:
    @inbounds ℓᶻ = - m.mixing_length.Cᴸᵏ * m.grid.zc[i]
    ℓᵇ = nan2inf(m.mixing_length.Cᴸᵇ * sqrt_e(m, i) / oncell(sqrt_∂B∂z, m, i))

    # Take hard minimum:
    ℓ = min(ℓᶻ, ℓᵇ)

    # Finally, limit by some factor of the local cell width
    ℓᵐⁱⁿ = m.mixing_length.Cᴸᵟ * Δf(m.grid, i)
    ℓ = max(ℓ, ℓᵐⁱⁿ)

    return ℓ
end

#
# Mixing length model from Tan et al 2018.
#

Base.@kwdef struct TanEtAl2018MixingLength{T} <: AbstractParameters
       Cᴸᵟ :: T = 1.0
       Cᴸᵏ :: T = 0.4
    Cᵃᵤₙₛ :: T = -100.0
    Cᵃₛₜₐ :: T = 2.7
    Cⁿᵤₙₛ :: T = 0.2   
    Cⁿₛₜₐ :: T = -1.0
end

@inline ζ(Ca, Qb, u★, z) = Ca * Qb / u★^3 * z
 
@inline function mixing_length(m::Model{<:TanEtAl2018MixingLength}, i)
    if isunstable(m)
        Cκ, Ca, Cn = m.mixing_length.Cᴸᵏ, m.mixing_length.Cᵃᵤₙₛ, m.mixing_length.Cⁿᵤₙₛ
        ℓᶻ = @inbounds Cκ * m.grid.zc[i] * (1 - ζ(Ca, m.state.Qb, u★(m)^3, m.grid.zc[i]))^Cn
    else
        Cκ, Ca, Cn = m.mixing_length.Cᴸᵏ, m.mixing_length.Cᵃₛₜₐ, m.mixing_length.Cⁿₛₜₐ
        ℓᶻ = @inbounds Cκ * m.grid.zc[i] * (1 - ζ(Ca, m.state.Qb, u★(m)^3, m.grid.zc[i]))^Cn
    end

    ℓᵟ = m.mixing_length.Cᴸᵟ * Δf(m.grid, i)
    ℓ = max(ℓᶻ, ℓᵟ)

    return ℓ
end

#
# Mixing length model due to Ignacio Lopez-Gomez + Clima
#

Base.@kwdef struct EquilibriumMixingLength{T} <: AbstractParameters
    Cᴸᵟ :: T = 1.0
    Cᴸᵏ :: T = 3.75
    Cᴸᵇ :: T = 0.64
end

"Returns τ = Cᴷ * ( (∂z U)² + (∂z V)^2 - Cᴾʳ * ∂z B ) at cell centers."
@inline tke_time_scale(m, i) = 
    m.tke_equation.Cᴷ * oncell(shear_squared, m, i) -
    m.tke_equation.Cᴷ * m.tke_equation.Cᴾʳᵩ * oncell_∂B∂z(m, i)

@inline function mixing_length(m::Model{<:EquilibriumMixingLength}, i)

    # Length scale associated with a steady TKE balance 
    τ = maxsqrt(tke_time_scale(m, i))
    ℓᵀᴷᴱ = sqrt_e(m, i) * sqrt(m.tke_equation.Cᴰ) / τ
    ℓᵀᴷᴱ = nan2inf(ℓᵀᴷᴱ)

    # Length scale associated with strongly-stratified turbulence
    ℓᵇ = m.mixing_length.Cᴸᵇ * sqrt_e(m, i) / maxsqrt(oncell_∂B∂z(m, i))
    ℓᵇ = nan2inf(ℓᵇ)

    # Length scale associated near-wall turbulence
    ℓʷ = @inbounds - m.mixing_length.Cᴸᵏ * m.grid.zc[i] * u★(m) / sqrt_e(m, i)

    # Hard minimum for now
    ℓ = min(ℓᵀᴷᴱ, ℓᵇ, ℓʷ)

    # Finally, limit by some factor of the local cell width
    ℓᵟ = m.mixing_length.Cᴸᵟ * Δf(m.grid, i)
    ℓ = max(ℓ, ℓᵟ)

    return ℓ
end
