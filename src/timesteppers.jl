# Equation abstraction and timesteppers for OceanTurb.jl

implicit_timestepping_methods = (:BackwardEuler,)

implicit(method) = method ∈ implicit_timestepping_methods

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

"Update the clock after one iteration."
function update!(clock, Δt)
  clock.time += Δt
  clock.iter += 1
  return nothing
end

#
# Forward Euler Timestepper
#

struct ForwardEulerTimestepper{T} <: Timestepper
  calc_rhs! :: Function
  rhs       :: T
  function ForwardEulerTimestepper(calc_rhs!, solution)
    rhs = deepcopy(solution)
    for fld in rhs
      set!(fld, 0)
    end
    new{typeof(rhs)}(calc_rhs!, rhs)
  end
end

function unpack(model::AbstractModel{TS}, i) where TS <: ForwardEulerTimestepper
  ϕ = model.solution[i]
  rhs = model.timestepper.rhs[i]
  return ϕ, rhs
end

# Forward Euler timestepping
function iterate!(model::AbstractModel{TS}, Δt) where TS <: ForwardEulerTimestepper

    # Evaluate the right-hand-side of ∂ϕ∂t for the current time-step.
    model.timestepper.calc_rhs!(model.timestepper.rhs, model)

    # Update the solution
    for j in eachindex(model.solution)
        @inbounds ϕ, rhs = unpack(model, j)
        for i in eachindex(ϕ)
            @inbounds ϕ[i] += Δt * rhs[i]
        end
    end

    update!(model.clock, Δt)

    return nothing
end

#
# Backward Euler Timestepper
#

#=
Could have

struct DiffusiveEquation{K}
    calc_rhs :: Function
    calc_lhs :: Function
    κ        :: K
end

and

function calc_lhs!(lhs, Δt, model::AbstractModel{E}) <: where E <: DiffusiveEquation
    calc_diffusive_lhs!(lhs, Δt, model.equation.κ, model)
    return nothing
end

But perhaps this is too complicated.

For now, our 'BackwardEulerTimestepper' can only treat diffusive terms
implicitly.
=#

struct BackwardEulerTimestepper{K, R, L} <: Timestepper
    calc_rhs! :: Function
    Κ         :: K
    rhs       :: R
    lhs       :: L
    function BackwardEulerTimestepper(calc_rhs!, Κ, solution, lhs)
        rhs = deepcopy(solution)
        for fld in rhs
            set!(fld, 0)
        end
        new{typeof(rhs), typeof(κ), typeof(lhs)}(calc_rhs!, Κ, rhs, lhs)
    end
end

function build_lhs_array(solution)
    lhs_array = []
    for fld in solution
        T = eltype(fld)
        A = arraytype(fld)
        N = length(fld)
        lhs_fld = Tridiagonal{T, A}(zeros(N-1), zeros(N), zeros(N-1))
        push!(lhs_array, lhs_fld)
    end
    return lhs_array
end

flux_div_op(m, K, face, cell) = K(m, face) / Δc(m, face) / Δf(m, cell)

# Build backward Euler operator for diffusive problems
function calc_diffusive_lhs!(lhs, Δt, K, m)
    @inbounds begin
        for j in eachindex(m.solution)
            ϕ = m.solution[j]
            Kϕ = K[j]

            for i in interior(ϕ)
                lhs.du[i] = -Δt*flux_div_op(m, Kϕ, i+1, i)
                lhs.d[i] = 1 + Δt*(flux_div_op(m, Kϕ, i+1, i) + flux_div_op(m, Kϕ, i, i))
                lhs.dl[i-1] = -Δt*flux_div_op(m, Kϕ, i, i)
            end

            # Bottom row
            be.lhs.du[1] = -Δt*flux_div_op(m, Kϕ, 2, 1)
            be.lhs.d[1] = 1 + Δt*flux_div_op(m, Kϕ, 2, 1)

            # Top row
            be.lhs.dl[length(ϕ)] = -Δt*flux_div_op(m, Kϕ, length(ϕ), length(ϕ))
            be.lhs.d[length(ϕ)]  = 1 + Δt*flux_div_op(m, Kϕ, length(ϕ), length(ϕ))
        end
    end

    return nothing
end

# Backward Euler timestepping for problems with diffusivity
function iterate!(m::AbstractModel{TS}, Δt) where TS <: BackwardEulerTimestepper

    # Evaluate the right-hand-side of ∂ϕ∂t for the current time-step.
    model.timestepper.calc_rhs!(model.timestepper.rhs, model)
    calc_diffusive_lhs!(model.timestepper.lhs, Δt, model.timestepper.K, model)

    # Update solution by inverting Tridiagonal lhs matrix
    for j in eachindex(m.solution)
        ϕdata = data(ϕ)
        @inbounds ϕdata .= data(rhs) \ be.lhs
    end

    update!(m.clock, Δt)

  return nothing
end
