Base.@kwdef struct TKEParameters{T} <: AbstractParameters
      Cᴰ :: T = 0.33  # Dissipation parameter
     Cᴷu :: T = 0.1   # Diffusivity parameter for velocity
    CᴷPr :: T = 0.74  # Diffusivity parameter for velocity
#     Cᴷe :: T = 0.1   # Ratio between turbulent kinetic energy and momentum diffusivity

    KU₀ :: T = 1e-6 # Interior viscosity for velocity
    KT₀ :: T = 1e-6 # Interior diffusivity for temperature
    KS₀ :: T = 1e-6 # Interior diffusivity for salinity
    Ke₀ :: T = 1e-6 # Interior diffusivity for salinity
end

# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.

@inline production(m, i) = KU(m, i) * oncell(shear_squared, m, i)

@inline ∂z_T(m, i) = ∂z(m.solution.T, i)
@inline ∂z_S(m, i) = ∂z(m.solution.S, i)

# Fallback for linear equations of state.
@inline buoyancy_flux(m, i) = - m.constants.g * (  m.constants.α * KT(m, i) * oncell(∂z_T, m, i)
                                                 - m.constants.β * KS(m, i) * oncell(∂z_S, m, i))

@inline dissipation(m, i) = @inbounds m.tke_equation.Cᴰ * sqrt_e(m, i)^3 / dissipation_length(m, i)
