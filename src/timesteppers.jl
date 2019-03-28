# Equation abstraction and timesteppers for OceanTurb.jl

import LinearAlgebra: Tridiagonal

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
    calc_rhs!   :: Function
    diffusivity :: K
    rhs         :: R
    lhs         :: L
    function BackwardEulerTimestepper(calc_rhs!, diffusivity, solution, lhs)
        rhs = deepcopy(solution)
        for fld in rhs
            set!(fld, 0)
        end
        new{typeof(diffusivity), typeof(rhs), typeof(lhs)}(calc_rhs!, diffusivity, rhs, lhs)
    end
end

function Tridiagonal(fld::Field)
    T = eltype(fld)
    A = arraytype(fld)
    N = length(fld)
    return Tridiagonal{T, A}(zeros(N-1), zeros(N), zeros(N-1))
end

function build_lhs(solution)
    lhs_1 = Tridiagonal(solution[1])
    lhs = [lhs_1]
    for i = 2:length(solution)
        lhs_i = Tridiagonal(solution[i])
        push!(lhs, lhs_i)
    end
    return lhs
end

flux_div_op(m, K::Function, face, cell) = K(m, face) / Δc(m.grid, face) / Δf(m.grid, cell)
flux_div_op(m, K::Number, face, cell) = K / Δc(m.grid, face) / Δf(m.grid, cell)

# Build backward Euler operator for diffusive problems
function calc_diffusive_lhs!(lhs, Δt, K, m)
    for j in eachindex(m.solution)
        @inbounds begin
            ϕ = m.solution[j]
            Kϕ = K[j]
            L = lhs[j]
        end

        for i in interiorindices(ϕ)
            @inbounds begin
                L.du[i] = -Δt*flux_div_op(m, Kϕ, i+1, i)
                L.d[i] = 1 + Δt*(flux_div_op(m, Kϕ, i+1, i) + flux_div_op(m, Kϕ, i, i))
                L.dl[i-1] = -Δt*flux_div_op(m, Kϕ, i, i)
            end
        end

        # Bottom row
        @inbounds L.du[1] = -Δt*flux_div_op(m, Kϕ, 2, 1)
        @inbounds L.d[1] = 1 + Δt*flux_div_op(m, Kϕ, 2, 1)

        # Top row
        @inbounds L.dl[end] = -Δt*flux_div_op(m, Kϕ, length(ϕ), length(ϕ))
        @inbounds L.d[end]  = 1 + Δt*flux_div_op(m, Kϕ, length(ϕ), length(ϕ))
    end

    return nothing
end

# Backward Euler timestepping for problems with diffusivity
function iterate!(model::AbstractModel{TS}, Δt) where TS <: BackwardEulerTimestepper

    # Evaluate the right-hand-side of ∂ϕ∂t for the current time-step.
    model.timestepper.calc_rhs!(model.timestepper.rhs, model)
    calc_diffusive_lhs!(model.timestepper.lhs, Δt, model.timestepper.diffusivity, model)

    # Update solution by inverting Tridiagonal lhs matrix
    for (j, ϕ) in enumerate(model.solution)
        @inbounds lhs = model.timestepper.lhs[j]
        @inbounds rhs = model.timestepper.rhs[j]

        for i in eachindex(rhs)
            @inbounds rhs[i] = ϕ[i] + Δt*rhs[i]
        end

        ldiv!(data(ϕ), lhs, data(rhs))
    end

    update!(model.clock, Δt)

  return nothing
end
