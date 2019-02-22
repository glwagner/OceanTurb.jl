# Model and Clock types for OceanTurb.jl

export
  Model,
  Clock

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
