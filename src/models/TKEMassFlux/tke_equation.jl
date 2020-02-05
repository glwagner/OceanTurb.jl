Base.@kwdef struct TKEParameters{T} <: AbstractParameters
     Cᴰ :: T = 0.305  # Dissipation parameter
     Cᴷ :: T = 1.0    # Diffusivity parameter 
   Cᴾʳᵩ :: T = 1.0    # Ratio between temperature and momentum diffusivity
   Cᴾʳₑ :: T = 0.1    # Ratio between turbulent kinetic energy and momentum diffusivity

    KU₀ :: T = 1e-6   # Interior viscosity for velocity
    KT₀ :: T = 1e-7   # Interior diffusivity for temperature
    KS₀ :: T = 1e-9   # Interior diffusivity for salinity
    Ke₀ :: T = 1e-6   # Interior diffusivity for salinity
end

# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.

@inline production(m, i) = KU(m, i) * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2)

@inline buoyancy_flux(m, i) = - m.constants.g * (
      m.constants.α * KT(m, i) * ∂z(m.solution.T, i)
    - m.constants.β * KS(m, i) * ∂z(m.solution.S, i)
    )

@inline dissipation(m, i) =
    @inbounds m.tke_equation.Cᴰ * maxsqrt(m.solution.e[i])^3 / dissipation_length(m, i)

#
# Turbulent kinetic energy wall model
#

# Fallbacks
TKEBoundaryConditions(T, wall_model) = DefaultBoundaryConditions(T)

Base.@kwdef struct PrescribedNearWallTKE{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
    Cʷw★ :: T = 0.2
    Cʷz★ :: T = 0.4
end

@inline update_near_wall_tke!(m::Model{L, H, <:PrescribedNearWallTKE}) where {L, H} =
    m.solution.e[m.grid.N] = m.tke_wall_model.Cʷu★ * u★(m)^2

Base.@kwdef struct SurfaceValueScaling{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
    Cʷw★ :: T = 0.2
    Cʷz★ :: T = 0.4
end

@inline (boundary_tke::SurfaceValueScaling)(model) = boundary_tke.Cʷu★ * u★(model)^2

TKEBoundaryConditions(T, wall_model::SurfaceValueScaling) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), ValueBoundaryCondition(wall_model))

Base.@kwdef struct SurfaceFluxScaling{T} <: AbstractParameters
    Cʷu★ :: T = 3.0
end

@inline (boundary_tke::SurfaceFluxScaling)(model) = - boundary_tke.Cʷu★ * u★(model)^3 # source...

TKEBoundaryConditions(T, wall_model::SurfaceFluxScaling) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(wall_model))
