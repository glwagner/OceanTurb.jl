Base.@kwdef struct TKEParameters{T} <: AbstractParameters
    Cᴰ :: T = 1.891 # Dissipation parameter
end

# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.

@inline production(m, i) = KU(m, i) * oncell(shear_squared, m, i)

@inline ∂z_T(m, i) = ∂z(m.solution.T, i)
@inline ∂z_S(m, i) = ∂z(m.solution.S, i)

@inline function buoyancy_flux(m, i) 
    return m.constants.g * (  m.constants.α * ( - KT(m, i) * oncell(∂z_T, m, i) + NLᵀ(m, i) )
                            - m.constants.β * ( - KS(m, i) * oncell(∂z_S, m, i) + NLˢ(m, i) ) )
end

@inline dissipation(m, i) = @inbounds m.tke_equation.Cᴰ * sqrt_e(m, i)^3 / dissipation_length(m, i)
