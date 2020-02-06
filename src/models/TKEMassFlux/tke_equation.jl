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

"Returns the eddy diffusivity at cell centers."
@inline function K(m, i) 
    if i <= 1 || i >= m.grid.N+1
        return zero(eltype(m.grid))
    else
        return @inbounds m.tke_equation.Cᴷ * diffusivity_mixing_length(m, i) * sqrt_e(m, i)
    end
end

@inline production(m, i) = KU(m, i) * oncell(shear_squared, m, i)

@inline ∂z_T(m, i) = ∂z(m.solution.T, i)
@inline ∂z_S(m, i) = ∂z(m.solution.S, i)

# Fallback for linear equations of state.
@inline buoyancy_flux(m, i) = - m.constants.g * (
        m.constants.α * KT(m, i) * oncell(∂z_T, m, i)
      - m.constants.β * KS(m, i) * oncell(∂z_S, m, i)
    )

@inline dissipation(m, i) =
    @inbounds m.tke_equation.Cᴰ * maxsqrt(m.solution.e[i])^3 / dissipation_length(m, i)

#
# Turbulent kinetic energy wall models
#

# Fallbacks
TKEBoundaryConditions(T, wall_model) = DefaultBoundaryConditions(T)

#
# Prescribes the "near wall" turbulent kinetic energy 
# at the upper grid cell at i = N.
#

Base.@kwdef struct PrescribedNearWallTKE{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
    Cʷw★ :: T = 0.2
    Cʷz★ :: T = 0.4
end

@inline update_near_wall_tke!(m::Model{L, H, <:PrescribedNearWallTKE}) where {L, H} =
    m.solution.e[m.grid.N] = m.tke_wall_model.Cʷu★ * u★(m)^2

#
# Prescribes the surface value of turbulent kinetic energy 
# via a ValueBoundaryCondition.
#

Base.@kwdef struct SurfaceValueScaling{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
    Cʷw★ :: T = 0.2
    Cʷz★ :: T = 0.4
end

@inline (boundary_tke::SurfaceValueScaling)(model) = boundary_tke.Cʷu★ * u★(model)^2

TKEBoundaryConditions(T, wall_model::SurfaceValueScaling) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), ValueBoundaryCondition(wall_model))

#
# Prescribes the surface flux of turbulent kinetic energy 
# as propotional to u★^3
#

Base.@kwdef struct SurfaceFluxScaling{T} <: AbstractParameters
    Cʷu★ :: T = 3.0
end

@inline (boundary_tke::SurfaceFluxScaling)(model) = - boundary_tke.Cʷu★ * u★(model)^3 # source...

TKEBoundaryConditions(T, wall_model::SurfaceFluxScaling) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(wall_model))

#
# Prescribes the surface flux of turbulent kinetic energy 
# as propotional to sqrt(e) u★^2, where sqrt(e) is the near-wall
# turbulent velocity at the upper grid cell i = N.
# This is a model for turbulent production across the surface
#

Base.@kwdef struct SurfaceProductionModel{T} <: AbstractParameters
    Cʷu★ :: T = 3.0
end

@inline (boundary_tke::SurfaceProductionModel)(m) = @inbounds - boundary_tke.Cʷu★ * sqrt_e(m, m.grid.N) * u★(m)^2 # source...

TKEBoundaryConditions(T, wall_model::SurfaceProductionModel) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(wall_model))
