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
# Explicit and implicit kernels
#

"Evaluate the right-hand-side of ∂ϕ∂t for an explicit time-stepping method."
function calc_explicit_rhs!(rhs, eqn, solution, m)

    # Function barrier for better performance and forced specialization
    # on Kϕ and Rϕ
    function kernel!(rhs, ϕ, K::KF, R::RF, m) where {KF, RF}
        for i in eachindex(ϕ)
            @inbounds rhs[i] = ∇K∇ϕ(K(m, i+1), K(m, i), ϕ, i) + R(m, i)
        end
    end

    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
        ϕ = solution[j]
        rhsϕ = rhs[j]
        Kϕ = eqn.K[j]
        Rϕ = eqn.R[j]

        kernel!(rhsϕ, ϕ, Kϕ, Rϕ, m)
    end
    return nothing
end

"Evaluate the right-hand-side of ∂ϕ∂t for an implicit time-stepping method."
function calc_implicit_rhs!(rhs, eqn, solution, m)

    N = m.grid.N

    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
        ϕ = solution[j]
        rhsϕ = rhs[j]
        Kϕ = eqn.K[j]
        Rϕ = eqn.R[j]

        for i in eachindex(rhsϕ)
            @inbounds rhsϕ[i] = Rϕ(m, i)
        end

        @inbounds rhsϕ[N] = ∇K∇ϕ(Kϕ(m, N+1), 0, ϕ, N) + Rϕ(m, N)
        @inbounds rhsϕ[1] = ∇K∇ϕ(0, Kϕ(m, 1), ϕ, 1) + Rϕ(m, 1)
    end
    return nothing
end

@inline flux_div_op(m, K::Function, face, cell) = K(m, face) / Δc(m.grid, face) / Δf(m.grid, cell)
@inline flux_div_op(m, K::Number, face, cell) = K / Δc(m.grid, face) / Δf(m.grid, cell)

"Build backward Euler operator for diffusive problems."
function calc_diffusive_lhs!(lhs, K, solution, Δt::T, m) where T
    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
        ϕ = solution[j]
        Kϕ = K[j]
        L = lhs[j]

        for i in interiorindices(ϕ)
            @inbounds begin
                L.du[i]   = -Δt * flux_div_op(m, Kϕ, i+1, i)
                L.d[i]    = one(T) + Δt * (flux_div_op(m, Kϕ, i+1, i) + flux_div_op(m, Kϕ, i, i))
                L.dl[i-1] = -Δt * flux_div_op(m, Kϕ, i, i)
            end
        end

        # Bottom row
        @inbounds L.du[1] = -Δt*flux_div_op(m, Kϕ, 2, 1)
        @inbounds L.d[1] = one(T) + Δt*flux_div_op(m, Kϕ, 2, 1)

        # Top row
        @inbounds L.dl[end] = -Δt*flux_div_op(m, Kϕ, length(ϕ), length(ϕ))
        @inbounds L.d[end]  = one(T) + Δt*flux_div_op(m, Kϕ, length(ϕ), length(ϕ))
    end

    return nothing
end

#
# Timestepping methods
#

function update!(bcs, eqn, solution, m)
    N = m.grid.N
    eqn.update!(m)
    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
        ϕ = solution[j]
        Kϕ = eqn.K[j]
        bcsϕ = bcs[j]
        fill_bottom_ghost_cell!(bcsϕ.bottom, ϕ, Kϕ(m, 1), m)
        fill_top_ghost_cell!(bcsϕ.top, ϕ, Kϕ(m, N+1), m, N)
    end
    return nothing
end

#
# ForwardEuler timestepper
#

"""
    ForwardEulerTimestepper(eqn, solution)

Construct a `FowardEulerTimestepper`.
"""
struct ForwardEulerTimestepper{E, R} <: Timestepper
    eqn :: E
    rhs :: R
    function ForwardEulerTimestepper(eqn, solution, args...)
        rhs = deepcopy(solution)
        for fld in rhs
            set!(fld, 0)
        end
        new{typeof(eqn), typeof(rhs)}(eqn, rhs)
    end
end


"Take one forward Euler step."
function forward_euler_step!(rhs, solution, Δt)
    for j in 1:length(solution)
        ϕ = solution[j]
        rhsϕ = rhs[j]

        for i in eachindex(ϕ)
            @inbounds ϕ[i] += Δt * rhsϕ[i]
        end
    end
    return nothing
end

# Forward Euler timestepping
function iterate!(m::AbstractModel{TS}, Δt) where TS <: ForwardEulerTimestepper
    update!(m.bcs, m.timestepper.eqn, m.solution, m)
    calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.solution, m)
    forward_euler_step!(m.timestepper.rhs, m.solution, Δt)
    tick!(m.clock, Δt)
    return nothing
end

#
# BackwardEuler timestepper
#

struct BackwardEulerTimestepper{E, R, L} <: Timestepper
    eqn :: E
    rhs :: R
    lhs :: L
    function BackwardEulerTimestepper(eqn, solution, lhs)
        rhs = deepcopy(solution)
        for fld in rhs
            set!(fld, 0)
        end
        new{typeof(eqn), typeof(rhs), typeof(lhs)}(eqn, rhs, lhs)
    end
end

function Tridiagonal(fld::AbstractField)
    T = eltype(fld)
    A = arraytype(fld)
    N = length(fld)
    return Tridiagonal{T, A}(zeros(N-1), zeros(N), zeros(N-1))
end

"""
    build_lhs(solution)

Build Tridiagonal matrices corresponding to the LHS
for each of `solution`'s equations.
"""
function build_lhs(solution)
    return ntuple(Val(length(solution))) do i
        Tridiagonal(solution[i])
    end
end

"Update solution by inverting Tridiagonal lhs matrix."
function backward_euler_step!(rhs, lhs, solution, Δt)
    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
        lhsϕ = lhs[j]
        rhsϕ = data(rhs[j])
        ϕ = data(solution[j])

        for i in eachindex(ϕ)
            @inbounds rhsϕ[i] = ϕ[i] + Δt*rhsϕ[i]
        end

        ldiv!(ϕ, lhsϕ, rhsϕ)
    end
    return nothing
end

"Step forward `m` by `Δt` with the backward Euler method."
function iterate!(m::AbstractModel{TS}, Δt) where TS <: BackwardEulerTimestepper
    update!(m.bcs, m.timestepper.eqn, m.solution, m)
    calc_implicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.solution, m)
    calc_diffusive_lhs!(m.timestepper.lhs, m.timestepper.eqn.K, m.solution, Δt, m)
    backward_euler_step!(m.timestepper.rhs, m.timestepper.lhs, m.solution, Δt)
    tick!(m.clock, Δt)
    return nothing
end
