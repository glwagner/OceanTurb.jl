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

function iterate!(model, Δt, Nt)
    for step = 1:Nt
        iterate!(model, Δt)
    end
    return nothing
end

"""
    iterate!(model; Δt, Nt)

Step `model` forward in time for `Nt` steps with step size Δt.
"""
iterate!(model; Δt, Nt) = iterate!(model, Δt, Nt)

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

function prepare!(bcs, K, solution, update_eqn!, N, m)
    update_eqn!(m)
    for j in eachindex(solution)
        @inbounds begin
            ϕ = solution[j]
            Kϕ = K[j]
            bcsϕ = bcs[j]
        end
        fill_bottom_ghost_cell!(bcsϕ.bottom, ϕ, Kϕ(m, 1), m)
        fill_top_ghost_cell!(bcsϕ.top, ϕ, Kϕ(m, N+1), m, N)
    end
    return nothing
end

@inline explicit_rhs_kernel(ϕ, K, R, m, i) = ∇K∇c(K(m, i+1), K(m, i), ϕ, i) + R(m, i)
@inline implicit_rhs_top(ϕ, K, R, m) = ∇K∇c(K(m, m.grid.N+1), 0, ϕ, m.grid.N) + R(m, m.grid.N)
@inline implicit_rhs_bottom(ϕ, K, R, m) = ∇K∇c(0, K(m, 1), ϕ, 1) + R(m, 1)
@propagate_inbounds implicit_rhs_kernel!(rhs, ϕ, R, m, i) = rhs[i] = R(m, i)

"Evaluate the right-hand-side of ∂ϕ∂t for the current time-step."
function calc_explicit_rhs!(rhs, eqn, solution, m)
    for j in eachindex(solution)

        @inbounds begin
            ϕ = solution[j]
            rhsϕ = rhs[j]
            Kϕ = eqn.K[j]
            Rϕ = eqn.R[j]
        end

        for i in eachindex(rhsϕ)
            @inbounds rhsϕ[i] = explicit_rhs_kernel(ϕ, Kϕ, Rϕ, m, i)
        end
    end
    return nothing
end

function calc_implicit_rhs!(rhs, eqn, solution, m)
    N = m.grid.N
    for j in eachindex(solution)

        @inbounds begin
            ϕ = solution[j]
            rhsϕ = rhs[j]
            Kϕ = eqn.K[j]
            Rϕ = eqn.R[j]
        end

        for i in interiorindices(rhsϕ)
            @inbounds implicit_rhs_kernel!(rhsϕ, ϕ, Rϕ, m, i)
        end

        @inbounds rhsϕ[N] = implicit_rhs_top(ϕ, Kϕ, Rϕ, m)
        @inbounds rhsϕ[1] = implicit_rhs_bottom(ϕ, Kϕ, Rϕ, m)
    end

    return nothing
end

function forward_euler_update!(rhs, solution, Δt)
    # Take one forward Euler step
    for j in eachindex(solution)

        @inbounds begin
            ϕ = solution[j]
            rhsϕ = rhs[j]
        end

        for i in eachindex(ϕ)
            @inbounds ϕ[i] += Δt * rhsϕ[i]
        end
    end
    return nothing
end

# Forward Euler timestepping
function iterate!(m::AbstractModel{TS}, Δt) where TS <: ForwardEulerTimestepper
    prepare!(m.bcs, m.timestepper.eqn.K, m.solution, m.timestepper.eqn.update!, m.grid.N, m)
    calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.solution, m)
    forward_euler_update!(m.timestepper.rhs, m.solution, Δt)
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

function Tridiagonal(fld::AbstractField)
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
function calc_diffusive_lhs!(Δt, lhs, K, solution, m)
    for j in eachindex(solution)
        @inbounds begin
            ϕ = solution[j]
            Kϕ = K[j]
            L = lhs[j]
        end

        for i in interiorindices(ϕ)
            @inbounds begin
                L.du[i]   = -Δt * flux_div_op(m, Kϕ, i+1, i)
                L.d[i]    = 1 + Δt * (flux_div_op(m, Kϕ, i+1, i) + flux_div_op(m, Kϕ, i, i))
                L.dl[i-1] = -Δt * flux_div_op(m, Kϕ, i, i)
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

"Update solution by inverting Tridiagonal lhs matrix."
function backward_euler_update!(rhs, lhs, solution, Δt)
    for j in eachindex(solution)
        @inbounds begin
            lhsϕ = lhs[j]
            rhsϕ = data(rhs[j])
            ϕ = data(solution[j])
        end

        for i in eachindex(rhsϕ)
            @inbounds rhsϕ[i] = ϕ[i] + Δt*rhsϕ[i]
        end

        ldiv!(ϕ, lhsϕ, rhsϕ)
    end
    return nothing
end

# Backward Euler timestepping for problems with diffusivity
function iterate!(m::AbstractModel{TS}, Δt) where TS <: BackwardEulerTimestepper
    prepare!(m.bcs, m.timestepper.eqn.K, m.solution, m.timestepper.eqn.update!, m.grid.N, m)
    calc_implicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.solution, m)
    calc_diffusive_lhs!(Δt, m.timestepper.lhs, m.timestepper.eqn.K, m.solution, m)
    backward_euler_update!(m.timestepper.rhs, m.timestepper.lhs, m.solution, Δt)
    tick!(m.clock, Δt)
    return nothing
end
