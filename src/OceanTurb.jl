module OceanTurb

export
  Model,
  Clock,
  Grid,
  Field,
  Constants,
  Parameters,
  Timestepper,

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

  # timesteppers.jl
  ForwardEulerTimestepper

using
  OffsetArrays

import Base: time

#
# Abstract types for OceanTurb.jl
#

abstract type Grid{T,A} end
abstract type Field{T,A,G} end
abstract type Constants end
abstract type Parameters end
abstract type Timestepper end

#
# Core OceanTurb.jl functionality
#

include("utils.jl")
include("grids.jl")
include("fields.jl")
include("timesteppers.jl")


mutable struct Clock{T}
  time::T
  iter::Int
end

Clock() = Clock(0.0, 0)

# Explain: what is a `Model`?

mutable struct Model{G<:Grid,C<:Constants,P<:Parameters,S,TS<:Timestepper,X,T<:AbstractFloat}
  grid::G
  constants::C
  parameters::P
  solution::S
  timestepper::TS
  auxiliaries::X # auxiliary variables needed to compute d(solution)/dt
  clock::Clock{T}
end

time(m::Model) = model.clock.time
iter(m::Model) = model.clock.iter

function Model(grid, constants, parameters, solution)
  Model(grid, constants, parameters, solution, NullTimestepper(), nothing, Clock())
end

function Model(grid, constants, parameters, solution, timestepper)
  Model(grid, constants, parameters, solution, timestepper, nothing, Clock())
end

#
# Ocean Turbulence Models
#

include("diffusion.jl")

end # module
