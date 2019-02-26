# Timesteppers for OceanTurb.jl

abstract type Timestepper end

"""
    Timestepper(stepper, args...)

Generalized timestepper constructor.
"""
function Timestepper(stepper, args...)
  fullsteppername = Symbol(stepper, :Timestepper)
  return eval(Expr(:call, fullsteppername, args...))
end

struct ForwardEulerTimestepper{S} <: Timestepper
  rhs::S
end

function step!(model, stepper::ForwardEulerTimestepper{S}, dt) where S <: Field
  calc_rhs!(stepper.rhs, equation, model)
  @. model.solution += dt*stepper.rhs
  nothing
end
