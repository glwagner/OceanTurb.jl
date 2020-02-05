Base.@kwdef struct SimpleMixingLength{T} <: AbstractParameters
    Cᶻ :: T = 0.4
    Cᵇ :: T = 0.1
    Cᵟ :: T = 1.0
end

@inline function diffusivity_mixing_length(m::Model{<:SimpleMixingLength}, i)
    # Two mixing lengths based on stratification and distance from surface:
    @inbounds ℓᶻ = - m.mixing_length.Cᶻ * m.grid.zf[i]
    ℓᵇ = nan2inf(m.mixing_length.Cᵇ * onface(sqrt_e, m, i) / maxsqrt(∂B∂z(m, i)))

    # Take hard minimum:
    ℓ = min(ℓᶻ, ℓᵇ)

    # Finally, limit by some factor of the local cell width
    ℓᵐⁱⁿ = m.mixing_length.Cᵟ * Δc(m.grid, i)
    ℓ = max(ℓ, ℓᵐⁱⁿ)

    return ℓ
end

@inline function dissipation_length(m::Model{<:SimpleMixingLength}, i)
    # Two mixing lengths based on stratification and distance from surface:
    @inbounds ℓᶻ = - m.mixing_length.Cᶻ * m.grid.zc[i]
    ℓᵇ = nan2inf(m.mixing_length.Cᵇ * sqrt_e(m, i) / oncell(sqrt_∂B∂z, m, i))

    # Take hard minimum:
    ℓ = min(ℓᶻ, ℓᵇ)

    # Finally, limit by some factor of the local cell width
    ℓᵐⁱⁿ = m.mixing_length.Cᵟ * Δf(m.grid, i)
    ℓ = max(ℓ, ℓᵐⁱⁿ)

    return ℓ
end

#
# Mixing length model from Tan et al 2018.
#

Base.@kwdef struct TanEtAl2018MixingLength{T} <: AbstractParameters
       Cᵟ :: T = 1.0
       Cᵏ :: T = 0.4
    Cᵃᵤₙₛ :: T = -100.0
    Cᵃₛₜₐ :: T = 2.7
    Cⁿᵤₙₛ :: T = 0.2   
    Cⁿₛₜₐ :: T = -1.0
end

@inline ζ(Ca, Qb, u★, z) = Ca * Qb / u★^3 * z

function diffusivity_mixing_length(m::Model{<:TanEtAl2018MixingLength}, i)

    if isunstable(m)
        Cκ, Ca, Cn = m.mixing_length.Cᵏ, m.mixing_length.Cᵃᵤₙₛ, m.mixing_length.Cⁿᵤₙₛ
        ℓᶻ = @inbounds Cκ * m.grid.zf[i] * (1 - ζ(Ca, m.state.Qb, u★(m)^3, m.grid.zf[i]))^Cn
    else
        Cκ, Ca, Cn = m.mixing_length.Cᵏ, m.mixing_length.Cᵃₛₜₐ, m.mixing_length.Cⁿₛₜₐ
        ℓᶻ = @inbounds Cκ * m.grid.zf[i] * (1 - ζ(Ca, m.state.Qb, u★(m)^3, m.grid.zf[i]))^Cn
    end

    ℓᵟ = m.mixing_length.Cᵟ * Δc(m.grid, i)
    ℓ = max(ℓᶻ, ℓᵟ)

    return ℓ
end

function dissipation_length(m::Model{<:TanEtAl2018MixingLength}, i)

    if isunstable(m)
        Cκ, Ca, Cn = m.mixing_length.Cᵏ, m.mixing_length.Cᵃᵤₙₛ, m.mixing_length.Cⁿᵤₙₛ
        ℓᶻ = @inbounds Cκ * m.grid.zc[i] * (1 - ζ(Ca, m.state.Qb, u★(m)^3, m.grid.zc[i]))^Cn
    else
        Cκ, Ca, Cn = m.mixing_length.Cᵏ, m.mixing_length.Cᵃₛₜₐ, m.mixing_length.Cⁿₛₜₐ
        ℓᶻ = @inbounds Cκ * m.grid.zc[i] * (1 - ζ(Ca, m.state.Qb, u★(m)^3, m.grid.zc[i]))^Cn
    end

    ℓᵟ = m.mixing_length.Cᵟ * Δf(m.grid, i)
    ℓ = max(ℓᶻ, ℓᵟ)

    return ℓ
end

#
# Mixing length model due to Ignacio Lopez + Clima
#

Base.@kwdef struct EquilibriumMixingLength{T} <: AbstractParameters
    Cᵟ :: T = 1.0
    Cᶻ :: T = 3.75
    Cᵇ :: T = 0.64
end

@inline tke_time_scale(m, i) = m.tke_equation.Cᴷᵤ * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2) -
                                m.tke_equation.Cᴷᵩ * ∂B∂z(m, i)

function diffusivity_mixing_length(m::Model{<:EquilibriumMixingLength}, i)

    # Length scale associated with a steady TKE balance 
    τ = maxsqrt(tke_time_scale(m, i))
    ℓᵀᴷᴱ = onface(sqrt_e, m, i) * sqrt(m.tke_equation.Cᴰ) / τ
    ℓᵀᴷᴱ = nan2inf(ℓᵀᴷᴱ)

    # Length scale associated with strongly-stratified turbulence
    ℓᵇ = m.mixing_length.Cᵇ * onface(sqrt_e, m, i) / maxsqrt(∂B∂z(m, i))
    ℓᵇ = nan2inf(ℓᵇ)

    # Length scale associated near-wall turbulence
    ℓʷ = @inbounds m.mixing_length.Cᶻ * m.grid.zf[i]

    # Hard minimum for now
    ℓ = min(ℓᵀᴷᴱ, ℓᵇ, ℓʷ)

    # Finally, limit by some factor of the local cell width
    ℓᵟ = m.mixing_length.Cᵟ * Δc(m.grid, i)
    ℓ = max(ℓ, ℓᵟ)

    return ℓ
end

function dissipation_length(m::Model{<:EquilibriumMixingLength}, i)

    # Length scale associated with a steady TKE balance 
    τ = maxsqrt(oncell(tke_time_scale, m, i))
    ℓᵀᴷᴱ = sqrt_e(m, i) * sqrt(m.tke_equation.Cᴰ) / τ
    ℓᵀᴷᴱ = nan2inf(ℓᵀᴷᴱ)

    # Length scale associated with strongly-stratified turbulence
    ℓᵇ = m.mixing_length.Cᵇ * sqrt_e(m, i) / oncell(sqrt_∂B∂z, m, i)
    ℓᵇ = nan2inf(ℓᵇ)

    # Length scale associated near-wall turbulence
    ℓʷ = @inbounds - m.mixing_length.Cᶻ * m.grid.zc[i] * u★(m) / sqrt_e(m, i)

    # Hard minimum for now
    ℓ = min(ℓᵀᴷᴱ, ℓᵇ, ℓʷ)

    # Finally, limit by some factor of the local cell width
    ℓᵟ = m.mixing_length.Cᵟ * Δf(m.grid, i)
    ℓ = max(ℓ, ℓᵟ)

    return ℓ
end
