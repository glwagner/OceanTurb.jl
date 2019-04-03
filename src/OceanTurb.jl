module OceanTurb

export # This file, core functionality:
    AbstractModel,
    AbstractSolution,
    AbstractParameters,
    Clock,
    Constants,
    time,
    iter,
    set!,
    reset!,

    # utils.jl
    second,
    minute,
    hour,
    day,
    year,
    stellaryear,
    Ω,
    pressenter,
    @zeros,
    @solution,
    @named_solution,
    @typed_solution,
    @pair_typed_solution,
    @add_standard_model_fields,
    @add_clock_grid_timestepper,

    # grids.jl
    Grid,
    UniformGrid,
    height,
    length,
    size,

    # fields.jl
    FieldLocation,
    Cell,
    Face,
    AbstractField,
    Field,
    CellField,
    FaceField,
    arraytype,
    data,
    nodes,
    Δc,
    Δf,
    onface,
    oncell,
    ∂z,
    ∂²z,
    ∂z!,
    set!,
    interiorindices,
    top,
    bottom,
    integral,
    top_flux_div,
    bottom_flux_div,

    # timesteppers.jl
    Equation,
    Timestepper,
    implicit,
    iterate!,
    ForwardEulerTimestepper,
    BackwardEulerTimestepper,

    # boundary_conditions.jl
    Flux,
    Gradient,
    Value,
    BoundaryCondition,
    FieldBoundaryConditions,
    FluxBoundaryCondition,
    ZeroFluxBoundaryConditions,
    ValueBoundaryCondition,
    GradientBoundaryCondition,
    set_top_bc!,
    set_bottom_bc!,
    set_top_flux_bc!,
    set_bottom_flux_bc!,
    set_flux_bcs!,
    set_bcs!,
    getbc,
    update_top_ghost_cell!,
    update_bottom_ghost_cell!,
    update_ghost_cells!,

    # Ocean turbulence models
    Diffusion,
    KPP,
    PacanowskiPhilander,
    ContinuousAdjustment,
    EddyDiffusivityMassFlux

using
    StaticArrays,
    OffsetArrays,
    LinearAlgebra

import Base: time, setproperty!

#
# Preliminary abstract types for OceanTurb.jl
#

abstract type AbstractParameters end
abstract type AbstractEquation end
abstract type Grid{T, A<:AbstractArray} end
abstract type Timestepper end
abstract type AbstractField{A<:AbstractArray, G<:Grid} end
abstract type AbstractSolution{N, T} <: FieldVector{N, T} end
abstract type AbstractModel{T, G, TS} end  # Explain: what is a `Model`?

#
# Core OceanTurb.jl functionality
#

include("utils.jl")
include("solvers.jl")
include("grids.jl")
include("boundary_conditions.jl")
include("fields.jl")
include("timesteppers.jl")

mutable struct Clock{T}
  time :: T
  iter :: Int
end

Clock() = Clock(0.0, 0)

"Get the current simulation time of the model."
time(m::AbstractModel) = m.clock.time

"Get the current iteration of the model."
iter(m::AbstractModel) = m.clock.iter

function reset!(clock)
  clock.time = 0
  clock.iter = 0
  return nothing
end

#
# Sugary things for solutions and models
#

"""
    set!(solution, kwargs...)

Set the fields of a solution. For example, use

T0 = rand(4)
S0(z) = exp(-z^2/10)
set!(solution, T=T0, S=S0)

To set solution.T and solution.S to T0 and S0.
"""
function set!(solution::AbstractSolution; kwargs...)
  for (k, v) in kwargs
    setproperty!(solution, k, v)
  end
  return nothing
end

function setproperty!(sol::AbstractSolution, c::Symbol, data::Union{Number, AbstractArray, Function})
  set!(getproperty(sol, c), data)
  return nothing
end

#
# Physical oceanic constants


struct Constants{T}
    g  :: T # Gravitiational acceleration
    cP :: T # Heat capacity of water
    ρ₀ :: T # Reference density
    α  :: T # Thermal expansion coefficient
    β  :: T # Haline expansion coefficient
    f  :: T # Coriolis parameter
end

function Constants(T=Float64; α=2.5e-4, β=8e-5, ρ₀=1035, cP=3992, f=0, g=9.81)
    Constants{T}(g, cP, ρ₀, α, β, f)
end

#
# Ocean Turbulence Models
#

include("models/diffusion.jl")
include("models/k_profile_parameterization.jl")
include("models/pacanowski_philander.jl")
include("models/continuous_convective_adjustment.jl")
include("models/eddy_diffusivity_mass_flux.jl")

end # module
