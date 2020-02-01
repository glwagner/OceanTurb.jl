Base.@kwdef struct TKEParameters{T} <: AbstractParameters

        CDe :: T = 2.0    # Dissipation parameter

       CK_U :: T = 0.1    # Diffusivity parameter
       CK_T :: T = 0.1    # Diffusivity parameter
       CK_e :: T = 0.1    # Diffusivity parameter

    Ca_unst :: T = -100.0
    Cb_unst :: T = -0.2

    Ca_stab :: T = 2.7
    Cb_stab :: T = -1.0

        KU₀ :: T = 1e-6   # Interior viscosity for velocity
        KT₀ :: T = 1e-7   # Interior diffusivity for temperature
        KS₀ :: T = 1e-9   # Interior diffusivity for salinity
        Ke₀ :: T = 1e-6   # Interior diffusivity for salinity
end

# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

"Return the turbuent velocity scale associated with convection."
@inline ωb(Qb, h) = abs(h * Qb)^(1/3)
@inline ωb(m::AbstractModel) = ωb(m.state.Qb, m.state.h)

@inline production(m, i) = KU(m, i) * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2)

@inline buoyancy_flux(m, i) = - m.constants.g * (
      m.constants.α * KT(m, i) * ∂z(m.solution.T, i)
    - m.constants.β * KS(m, i) * ∂z(m.solution.S, i)
    )

@inline dissipation(m, i) =
    @inbounds m.tke.CDe * maxsqrt(m.solution.e[i])^3 / mixing_length(m, i)

@inline maxsqrt(ϕ::T) where T = sqrt(max(zero(T), ϕ))
@inline maxsqrt(ϕ, i) = @inbounds sqrt(max(0, ϕ[i]))
