# Fallbacks
@inline diffusivity_mixing_length(m, i) = mixing_length(m, i)
@inline dissipation_length(m, i) = mixing_length(m, i)

#
# A "simple" mixing length
#

Base.@kwdef struct SimpleMixingLength{T} <: AbstractParameters
    Cᴸᵇ :: T = 3.9302
end

@inline function mixing_length(m::Model{<:SimpleMixingLength}, i)
    # Two mixing lengths based on stratification and distance from surface:
    @inbounds ℓᶻ = - m.grid.zc[i]
    ℓᵇ = nan2inf(m.mixing_length.Cᴸᵇ * sqrt_e(m, i) / oncell(sqrt_∂B∂z, m, i))

    # Take hard minimum:
    ℓ = min(ℓᶻ, ℓᵇ)

    return ℓ
end

#
# Mixing length model due to Ignacio Lopez-Gomez + Clima
#

Base.@kwdef struct EquilibriumMixingLength{T} <: AbstractParameters
    Cᴸᵟ :: T = 1.0
    Cᴸʷ :: T = 1.688 # Limits to Von-Karman constant for stress-driven turbulence
    Cᴸᵇ :: T = 1.664
end

@inline function mixing_length(m::Model{<:EquilibriumMixingLength}, i)
    @inbounds e = m.solution.e[i] # TKE at cell i

    # "Smallest" diffusivity limited by grid spacing.
    ℓᵟ = m.mixing_length.Cᴸᵟ * Δf(m.grid, i)

    if e <= 0 # shortcut when TKE is zero (diffusivity is zero in this case regardless).
        return ℓᵟ
    else
        # For notational convenience
        Cᴸᵇ = m.mixing_length.Cᴸᵇ
        Cᴸʷ = m.mixing_length.Cᴸʷ
        Cᴰ  = m.tke_equation.Cᴰ

        # Extract buoyancy frequency and shear squared.
        N² = oncell_∂B∂z(m, i) 
        S² = oncell(shear_squared, m, i)

        # Length scale associated with a production-buoyancy flux-dissipation balance
        # in the TKE budget. Valid only for a linear equation of state with a shared diffusivity
        # for salinity and temperature.
        #
        # Note that ω² = Cᴷu * S² - Cᴷc * N² can be negative, in which case this model predicts ℓᵀᴷᴱ=0.
        ω² = Cᴷu(m, i) * S² - Cᴷc(m, i) * N²
        ℓᵀᴷᴱ = maxsqrt(Cᴰ * e / ω²)

        # Length-scale limitation by strong stratification.
        ℓᵇ = N² <= 0 ? Inf : Cᴸᵇ * √(e / N²)

        # Near-wall length scale:
        ℓʷ = @inbounds - Cᴸʷ * m.grid.zc[i] * u★(m) / √e

        # Hard minimum for now
        ℓ = min(ℓᵀᴷᴱ, ℓᵇ, ℓʷ)

        return max(ℓ, ℓᵟ) # limits mixing length to be larger than grid spacing
    end
end
