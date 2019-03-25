module OceanTurb

export # This file, core functionality:
    AbstractModel,
    AbstractSolution,
    AbstractParameters,
    Clock,
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
    @specify_solution,
    @add_standard_model_fields,

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
    interior,
    top,
    bottom,

    # equations.jl
    Equation,

    # timesteppers.jl
    Timestepper,
    iterate!,
    unpack,
    ForwardEulerTimestepper,

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

    # Ocean turbulence models
    Diffusion,
    KPP,
    PacanowskiPhilander

using
    StaticArrays,
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
abstract type AbstractModel{TS<:Timestepper, G<:Grid, E<:AbstractEquation, T<:AbstractFloat} end  # Explain: what is a `Model`?

#
# Core OceanTurb.jl functionality
#

include("utils.jl")
include("grids.jl")
include("fields.jl")
include("boundary_conditions.jl")
include("operators.jl")
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
#

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

include("diffusion.jl")
include("k_profile_parameterization.jl")
include("pacanowski_philander.jl")

end # module
