module OceanTurb

export
  Model,
  Clock,
  Parameters

using
  OffsetArrays

#
# Core OceanTurb.jl functionality
#

include("utils.jl")
include("grids.jl")
include("fields.jl")
include("equations.jl")
include("timesteppers.jl")

abstract type Parameters end

mutable struct Clock{T}
  t::T
  step::Int
end

Clock() = Clock(0.0, 0)

# Explain: what is a `Model`?

mutable struct Model{G,C,P,E,S,X,T}
  grid::G
  constants::C
  parameters::P
  equation::E
  solution::S
  auxiliaries::X # auxiliary variables
  clock::Clock{T}
end

function Model(grid, constants, parameters, equation, solution)
  Model(grid, constants, parameters, equation, solution, nothing, Clock())
end


#
# Ocean Turbulence Models
#

include("models/diffusion.jl")

end # module
