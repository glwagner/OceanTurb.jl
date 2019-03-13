# Timesteppers for OceanTurb.jl

"""
    Timestepper(stepper, args...)

Generalized timestepper constructor.
"""
function Timestepper(stepper, args...)
  fullsteppername = Symbol(stepper, :Timestepper)
  return eval(Expr(:call, fullsteppername, args...))
end

"""
    iterate!(model, Δt, nt=1)

Step `model` forward in time by one time-step with step-size `Δt`.
"""
function iterate!(model, Δt, nt)
  for step = 1:nt
    iterate!(model, Δt)
  end
  return nothing
end

function unpack(model, i)
  ϕ = model.solution[i]
  ∂ϕ∂t = model.equation[i]
  rhs = model.timestepper.rhs[i]
  bcs = model.bcs[i]
  return ϕ, ∂ϕ∂t, rhs, bcs
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

function update!(clock, Δt)
  clock.time += Δt
  clock.iter += 1
  return nothing
end

# Forward Euler timestepping
function iterate!(model::AbstractModel{TS}, Δt) where TS <: ForwardEulerTimestepper

  for j in eachindex(model.solution)
    ϕ, ∂ϕ∂t, rhs, bcs = unpack(model, j)

    # Interior step
    for i in interior(ϕ)
      @inbounds rhs.data[i] = ∂ϕ∂t(model, i)
    end

    # Boundary conditions
    rhs.data[end] = ∂ϕ∂t(model, bcs.top)
    rhs.data[1]   = ∂ϕ∂t(model, bcs.bottom)
  end

  # Update solution
  for j in eachindex(model.solution)
    ϕ, ∂ϕ∂t, rhs, bcs = unpack(model, j)
    @. ϕ.data += Δt*rhs.data
  end

  update!(model.clock, Δt)

  return nothing
end
