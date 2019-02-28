# Timesteppers for OceanTurb.jl

"""
    Timestepper(stepper, args...)

Generalized timestepper constructor.
"""
function Timestepper(stepper, args...)
  fullsteppername = Symbol(stepper, :Timestepper)
  return eval(Expr(:call, fullsteppername, args...))
end

struct NullTimestepper <: Timestepper; end

struct ForwardEulerTimestepper{S} <: Timestepper
  rhs::S
  function ForwardEulerTimestepper(solution)
    rhs = similar(solution)
    return new{typeof(rhs)}(rhs)
  end
end


function step!(model, stepper::ForwardEulerTimestepper{S}, dt) where S <: Field
  calc_rhs!(stepper.rhs, equation, model)
  @. model.solution += dt*stepper.rhs
  nothing
end
