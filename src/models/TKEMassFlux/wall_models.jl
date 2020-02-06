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

Base.@kwdef struct SurfaceTKEProductionModel{T} <: AbstractParameters
    Cʷu★ :: T = 3.0
end

@inline (boundary_tke::SurfaceTKEProductionModel)(model) = - boundary_tke.Cʷu★ * u★(model)^3 # source...

TKEBoundaryConditions(T, wall_model::SurfaceTKEProductionModel) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(wall_model))

#
# Prescribes the surface flux of turbulent kinetic energy 
# as propotional to sqrt(e) u★^2, where sqrt(e) is the near-wall
# turbulent velocity at the upper grid cell i = N.
#

Base.@kwdef struct MixedSurfaceTKEProductionModel{T} <: AbstractParameters
    Cʷu★ :: T = 1.0
end

@inline (boundary_tke::MixedSurfaceTKEProductionModel)(m) = @inbounds - boundary_tke.Cʷu★ * sqrt_e(m, m.grid.N) * u★(m)^2 # source...

TKEBoundaryConditions(T, wall_model::MixedSurfaceTKEProductionModel) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(wall_model))
