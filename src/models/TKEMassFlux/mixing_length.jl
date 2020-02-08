# Fallbacks
@inline diffusivity_mixing_length(m, i) = mixing_length(m, i)
@inline dissipation_length(m, i) = mixing_length(m, i)

"Returns (∂z U)^2 + (∂z V)^2 at cell interfaces."
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
# Mixing length model due to Ignacio Lopez-Gomez + Clima
#

Base.@kwdef struct EquilibriumMixingLength{T} <: AbstractParameters
    Cᴸᵟ :: T = 1.0
    Cᴸᵏ :: T = 0.4
    Cᴸᵇ :: T = 0.5
end

"Returns τ² = 1 / Cᴷ * ( (∂z U)² + (∂z V)^2 - Cᴾʳ * ∂z B ) at cell centers."
@inline function tke_time_scale_squared(m, i)
    Cᴷᵤ, Cᴾʳ = m.tke_equation.Cᴷᵤ, m.tke_equation.Cᴾʳ
    ω² = Cᴷᵤ * (oncell(shear_squared, m, i) - Cᴾʳ * oncell_∂B∂z(m, i))
    return 1 / ω²
end

@inline function mixing_length(m::Model{<:EquilibriumMixingLength}, i)

    @inbounds e = m.solution.e[i] # TKE at cell i

    # "Smallest" diffusivity limited by grid spacing.
    ℓᵟ = m.mixing_length.Cᴸᵟ * Δf(m.grid, i)

    if e <= 0 # shortcut when TKE is zero (diffusivity is zero in this case regardless).
        return ℓᵟ
    else
        # For notational convenience
        Cᴸᵇ, Cᴸᵏ, = m.mixing_length.Cᴸᵟ, m.mixing_length.Cᴸᵇ, m.mixing_length.Cᴸᵏ
        Cᴰ = m.tke_equation.Cᴰ

        # Length scale associated with a production-buoyancy flux-dissipation balance
        # in the TKE budget.
        τ = tke_time_scale(m, i)
        ℓᵀᴷᴱ = √(Cᴰ * τ² * e)

        # Length-scale limitation by strong stratification.
        N = oncell_∂B∂z(m, i) 
        ℓᵇ = N == 0 ? Inf : Cᴸᵇ * √e / N

        # Near-wall length scale:
        ℓʷ = @inbounds - Cᴸᵏ * m.grid.zc[i] * u★(m) / √e

        # Hard minimum for now
        ℓ = min(ℓᵀᴷᴱ, ℓᵇ, ℓʷ)

        return max(ℓ, ℓᵟ) # limits mixing length to be larger than grid spacing
    end
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

