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

struct ForwardEulerTimestepper{T} <: Timestepper
  rhs::T
  function ForwardEulerTimestepper(solution)
    rhs = deepcopy(solution)
    for fld in rhs
      set!(fld, 0)
    end
    new{typeof(rhs)}(rhs)
  end
end
