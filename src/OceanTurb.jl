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
    @use_pyplot_utils,

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
    @prefixed_solution,
    @typed_solution,
    @pair_typed_solution,
    @add_standard_model_fields,
    @add_clock_grid_timestepper,
    run_until!,

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
    absolute_error,
    relative_error,

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
    BoundaryConditions,
    FluxBoundaryCondition,
    ValueBoundaryCondition,
    GradientBoundaryCondition,
    FieldBoundaryConditions,
    ZeroFluxBoundaryConditions,
    DefaultBoundaryConditions,
    set_bcs!,
    getbc,
    fill_top_ghost_cell!,
    fill_bottom_ghost_cell!,
    fill_ghost_cells!,

    # Ocean turbulence models
    Diffusion,
    KPP,
    ModularKPP,
    PacanowskiPhilander,
    KPP_TKE

using
    Statistics,
    StaticArrays,
    OffsetArrays,
    LinearAlgebra

import Base: time, setproperty!

#
# Preliminary abstract types for OceanTurb.jl
#

abstract type AbstractParameters end
abstract type AbstractEquation end
abstract type Grid{T, A} end
abstract type Timestepper end
abstract type AbstractField{A, G, T} end
abstract type AbstractSolution{N, T} <: FieldVector{N, T} end
abstract type AbstractModel{T, G, TS} end

#
# Core OceanTurb.jl functionality
#

include("utils.jl")
include("solvers.jl")
include("grids.jl")
include("boundary_conditions.jl")
include("fields.jl")
include("equations.jl")
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

"""
    Constants(T=Float64; α=2.5e-4, β=8e-5, ρ₀=1035, cP=3992, f=0, g=9.81)

Construct `Constants` with
    * thermal expansion coefficient `α` [C⁻¹]
    * haline contraction coefficient `β` [psu⁻¹]
    * reference density `ρ₀` [kg m⁻³]
    * heat capacity `cP` [...]
    * Coriolis parameter `f` [s⁻¹]
    * gravitational acceleration `g` [m² s⁻¹]
"""
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

include("models/Diffusion.jl")
include("models/KPP/KPP.jl")
include("models/ModularKPP/ModularKPP.jl")
include("models/PacanowskiPhilander.jl")

# Convenient utilities for plotting
macro use_pyplot_utils()
    return esc(quote
        using PyPlot, PyCall
        include(joinpath(@__DIR__, "..", "plotting", "pyplot_utils.jl"))
        using Main.OceanTurbPyPlotUtils
    end
    )
end

end # module
