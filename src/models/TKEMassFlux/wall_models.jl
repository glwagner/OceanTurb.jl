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
end

@inline update_near_wall_tke!(m::Model{L, H, <:PrescribedNearWallTKE}) where {L, H} =
    m.solution.e[m.grid.N] = m.tke_wall_model.Cʷu★ * u★(m)^2

#
# Prescribes the surface value of turbulent kinetic energy 
# via a ValueBoundaryCondition.
#

Base.@kwdef struct PrescribedSurfaceTKE{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
end

@inline (boundary_tke::PrescribedSurfaceTKE)(model) = boundary_tke.Cʷu★ * u★(model)^2

TKESurfaceConditions(T, wall_model::PrescribedSurfaceTKE) =
    FieldSurfaceConditions(GradientSurfaceCondition(-zero(T)), ValueSurfaceCondition(wall_model))

#
# Prescribes the surface flux of turbulent kinetic energy 
# as propotional to u★^3
#

Base.@kwdef struct PrescribedSurfaceTKEFlux{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
end

@inline (boundary_tke::PrescribedSurfaceTKEFlux)(model) = - boundary_tke.Cʷu★ * u★(model)^3 # source...

TKEBoundaryConditions(T, wall_model::PrescribedSurfaceTKEFlux) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(wall_model))

#
# Prescribes the surface flux of turbulent kinetic energy 
# as propotional to sqrt(e) u★^2, where sqrt(e) is the near-wall
# turbulent velocity at the upper grid cell i = N.
#

Base.@kwdef struct DynamicSurfaceTKEFlux{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
end

@inline (boundary_tke::DynamicSurfaceTKEFlux)(m) = @inbounds - boundary_tke.Cʷu★ * sqrt_e(m, m.grid.N) * u★(m)^2 # source...

TKEBoundaryConditions(T, wall_model::DynamicSurfaceTKEFlux) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(wall_model))
