module OceanTurb

export # by this file:
  Model,
  Clock,
  Grid,
  Field,
  Constants,
  AbstractParameters,
  AbstractSolution,
  Timestepper,
  BoundaryCondition,
  AbstractModel,
  time,
  iter,
  set!,

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
  @def_solution_fields,

  # grids.jl
  UniformGrid,

  # fields.jl
  Field,
  CellField,
  FaceField,
  dzc,
  dzf,
  ∂z,
  ∂²z,
  ∂z!,
  set!,
  interior,
  top,
  bottom,

  # timesteppers.jl
  NullTimestepper,
  ForwardEulerTimestepper,

  # boundary_conditions.jl
  FieldBoundaryConditions,
  FluxBC,
  ConstBC

using
  StaticArrays,
  Parameters

import Base: time

#
# Abstract types for OceanTurb.jl
#

abstract type Grid{T,A} end
abstract type Field{A,G} end
abstract type AbstractSolution{N,T<:Field} <: FieldVector{N,T} end
abstract type AbstractParameters end
abstract type Constants end
abstract type Timestepper end
abstract type BoundaryCondition{side} end
abstract type AbstractModel{G<:Grid,TS<:Timestepper} end  # Explain: what is a `Model`?

#
# Core OceanTurb.jl functionality
#

include("utils.jl")
include("grids.jl")
include("fields.jl")
include("boundary_conditions.jl")
include("timesteppers.jl")

mutable struct Clock{T}
  time::T
  iter::Int
end

Clock() = Clock(0.0, 0)
time(m::AbstractModel) = model.clock.time
iter(m::AbstractModel) = model.clock.iter

#=
# Future functionality!
function set!(solution::FieldVector; kwargs...)
  for (k, v) in kwargs
    setproperty!(solution, k, v)
  end
  return nothing
end
=#

#
# Ocean Turbulence Models
#

include("diffusion.jl")

end # module
