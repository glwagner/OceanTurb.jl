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

Base.@kwdef struct PrescribedSurfaceTKEValue{T} <: AbstractParameters
    Cʷu★ :: T = 3.75
end

@inline (boundary_tke::PrescribedSurfaceTKEValue)(model) = boundary_tke.Cʷu★ * u★(model)^2

TKEValueSurfaceConditions(T, wall_model::PrescribedSurfaceTKEValue) =
    FieldSurfaceConditions(GradientSurfaceCondition(-zero(T)), ValueSurfaceCondition(wall_model))

#
# Prescribes the surface flux of turbulent kinetic energy 
# as propotional to the sum of u★^3 and wΔ^3 = Δz * max(0, Qb), where Qb is the buoyancy flux.
# The free parameters for both u★ and wΔ are multiplied by the TKE dissipation parameter
# to prevent their mutual correlation.
#

Base.@kwdef struct PrescribedSurfaceTKEFlux{T} <: AbstractParameters
    Cʷu★ :: T = 1.3717
    CʷwΔ :: T = 1.0
end

@inline prescribed_surface_tke_flux(model) =  - model.tke_equation.Cᴰ * (   model.tke_wall_model.Cʷu★ * u★(model)^3
                                                                          + model.tke_wall_model.CʷwΔ * wΔ³(model) )  

TKEBoundaryConditions(T, wall_model::PrescribedSurfaceTKEFlux) =
    FieldBoundaryConditions(GradientBoundaryCondition(-zero(T)), FluxBoundaryCondition(prescribed_surface_tke_flux))
