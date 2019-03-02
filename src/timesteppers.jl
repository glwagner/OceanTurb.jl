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

function step!(model, Δt, nt)
  for step = 1:nt
    step!(model, Δt)
  end
  return nothing
end

function unpack(model, i)
  c = model.solution[i]
  ∂c∂t = model.equation[i]
  rhs = model.timestepper.rhs[i]
  bcs = model.bcs[i]
  return c, ∂c∂t, rhs, bcs
end

#
# Forward Euler Timestepper
#

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

# Forward Euler timestepping
function step!(model, Δt)

  for i in eachindex(model.solution)
    c, ∂c∂t, rhs, bcs = unpack(model, i)

    # Interior step
    for i in interior(c)
      @inbounds rhs.data[i] = ∂c∂t(model, i)
    end

    # Boundary conditions
    rhs.data[end] = ∂c∂t(model, bcs.top)
    rhs.data[1]   = ∂c∂t(model, bcs.bottom)
  end

  # Time step
  for i in eachindex(model.solution)
    c, ∂c∂t, rhs, bcs = unpack(model, i)
    @. c.data += Δt*rhs.data
  end

  return nothing
end
