# Equation abstraction and timesteppers for OceanTurb.jl

import LinearAlgebra: Tridiagonal
import Base: @propagate_inbounds

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
function tick!(clock, Δt)
    clock.time += Δt
    clock.iter += 1
    return nothing
end

#
# ForwardEuler timestepper
#

struct ForwardEulerTimestepper{R, ER, EK} <: Timestepper
    eqn :: Equation{ER, EK}
    rhs :: R
    function ForwardEulerTimestepper(eqn::Equation{R, K}, solution, args...) where {R, K}
        rhs = deepcopy(solution)
        for fld in rhs
            set!(fld, 0)
        end
        new{typeof(rhs), R, K}(eqn, rhs)
    end
end

function update!(m)
    m.timestepper.eqn.update!(m)
    for j in eachindex(m.solution)
        @inbounds ϕ, rhsϕ, Rϕ, Kϕ, bcsϕ = unpack(m, j)
        @inbounds update_ghost_cells!(ϕ, Kϕ(m, 1), Kϕ(m, m.grid.N+1), m, bcsϕ)
    end
    return nothing
end

explicit_rhs_kernel(ϕ, Kϕ, Rϕ, m, i) = ∇K∇c(Kϕ(m, i+1), Kϕ(m, i), ϕ, i) + Rϕ(m, i)
explicit_rhs_kernel(ϕ, Kϕ, ::Nothing, m, i) = ∇K∇c(Kϕ(m, i+1), Kϕ(m, i), ϕ, i)

implicit_rhs_top(ϕ, Kϕ, Rϕ, m) = @inbounds ∇K∇c(Kϕ(m, m.grid.N+1), 0, ϕ, m.grid.N) + Rϕ(m, m.grid.N)
implicit_rhs_top(ϕ, Kϕ, ::Nothing, m) = @inbounds ∇K∇c(Kϕ(m, m.grid.N+1), 0, ϕ, m.grid.N)

implicit_rhs_bottom(ϕ, Kϕ, Rϕ, m) = @inbounds ∇K∇c(0, Kϕ(m, 1), ϕ, 1) + Rϕ(m, 1)
implicit_rhs_bottom(ϕ, Kϕ, ::Nothing, m) = @inbounds ∇K∇c(0, Kϕ(m, 1), ϕ, 1)

implicit_rhs_kernel!(rhs, ϕ, R, m, i) = rhs[i] = R(m, i)
implicit_rhs_kernel!(rhs, ϕ, ::Nothing, m, i) = nothing

"Evaluate the right-hand-side of ∂ϕ∂t for the current time-step."
function calc_explicit_rhs!(rhs, eqn, m)
    for (j, rhsϕ) in enumerate(rhs)
        @inbounds ϕ, rhsϕ, Rϕ, Kϕ, bcsϕ = unpack(m, j)
        for i in eachindex(rhsϕ)
            @inbounds rhsϕ[i] = explicit_rhs_kernel(ϕ, Kϕ, Rϕ, m, i)
        end
    end
    return nothing
end

function calc_implicit_rhs!(rhs, eqn, m)
    N = m.grid.N
    for (j, rhsϕ) in enumerate(rhs)
        @inbounds ϕ, rhsϕ, Rϕ, Kϕ, bcsϕ = unpack(m, j)

        for i in interiorindices(rhsϕ)
            @inbounds implicit_rhs_kernel!(rhsϕ, ϕ, Rϕ, m, i)
        end

        @inbounds rhsϕ[N] = implicit_rhs_top(ϕ, Kϕ, Rϕ, m)
        @inbounds rhsϕ[1] = implicit_rhs_bottom(ϕ, Kϕ, Rϕ, m)
    end

    return nothing
end

# Forward Euler timestepping
function iterate!(m::AbstractModel{TS}, Δt) where TS <: ForwardEulerTimestepper

    update!(m)
    calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m)

    # Take one forward Euler step
    for j in eachindex(m.solution)
        @inbounds ϕ, rhsϕ, Rϕ, Kϕ, bcsϕ = unpack(m, j)
        for i in eachindex(ϕ)
            @inbounds ϕ[i] += Δt * rhsϕ[i]
        end
    end

    tick!(m.clock, Δt)

    return nothing
end


#
# BackwardEuler timestepper
#

struct BackwardEulerTimestepper{R, L, ER, EK} <: Timestepper
    eqn :: Equation{ER, EK}
    rhs :: R
    lhs :: L
    function BackwardEulerTimestepper(eqn::Equation{R, K}, solution, lhs) where {R, K}
        rhs = deepcopy(solution)
        for fld in rhs
            set!(fld, 0)
        end
        new{typeof(rhs), typeof(lhs), R, K}(eqn, rhs, lhs)
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
function iterate!(m::AbstractModel{TS}, Δt) where TS <: BackwardEulerTimestepper

    # Evaluate the right-hand-side of ∂ϕ∂t for the current time-step.
    update!(m)
    calc_implicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m)
    calc_diffusive_lhs!(m.timestepper.lhs, Δt, m.timestepper.eqn.K, m)

    # Update solution by inverting Tridiagonal lhs matrix
    for (j, ϕ) in enumerate(m.solution)
        @inbounds lhs = m.timestepper.lhs[j]
        @inbounds rhs = m.timestepper.rhs[j]

        for i in eachindex(rhs)
            @inbounds rhs[i] = ϕ[i] + Δt*rhs[i]
        end

        ldiv!(data(ϕ), lhs, data(rhs))
    end

    tick!(m.clock, Δt)

  return nothing
end

function unpack(model::AbstractModel{TS},
    i) where TS <: Union{ForwardEulerTimestepper, BackwardEulerTimestepper}

      ϕ = model.solution[i]
    rhs = model.timestepper.rhs[i]
      R = model.timestepper.eqn.R[i]
      K = model.timestepper.eqn.K[i]
    bcs = model.bcs[i]

    return ϕ, rhs, R, K, bcs
end
