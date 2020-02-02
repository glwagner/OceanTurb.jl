Base.@kwdef struct TKEParameters{T} <: AbstractParameters
        CDe :: T = 0.305  # Dissipation parameter

       CK_U :: T = 1.0    # Diffusivity parameter
       CK_T :: T = 2.0    # Diffusivity parameter
       CK_e :: T = 0.1    # Diffusivity parameter

        KU₀ :: T = 1e-6   # Interior viscosity for velocity
        KT₀ :: T = 1e-7   # Interior diffusivity for temperature
        KS₀ :: T = 1e-9   # Interior diffusivity for salinity
        Ke₀ :: T = 1e-6   # Interior diffusivity for salinity
end

# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

@inline production(m, i) = KU(m, i) * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2)

@inline buoyancy_flux(m, i) = - m.constants.g * (
      m.constants.α * KT(m, i) * ∂z(m.solution.T, i)
    - m.constants.β * KS(m, i) * ∂z(m.solution.S, i)
    )

@inline dissipation(m, i) =
    @inbounds m.tke_equation.CDe * maxsqrt(m.solution.e[i])^3 / mixing_length_cell(m, i)

