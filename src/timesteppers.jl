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

function time_step!(model, Δt, Nt)
    for step = 1:Nt
        time_step!(model, Δt)
    end
    return nothing
end

"""
    time_step!(model; Δt, Nt)

Step `model` forward in time for `Nt` steps with step size Δt.
"""
time_step!(model; Δt, Nt) = time_step!(model, Δt, Nt)

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
function calc_explicit_rhs!(rhs, eqn, bcs, solution, m)

    # Function barrier for better performance and forced specialization
    # on Kϕ and Rϕ
    function kernel!(rhs, bcs, ϕ, A::AF, L::LF, K::KF, R::RF, m) where {AF, LF, KF, RF}
        for i in eachindex(ϕ)
            @inbounds rhs[i] = (-L(m, i) * ϕ[i] - ∂zA(A(m, i+1), A(m, i), ϕ, i)
                                    + ∇K∇ϕ(K(m, i+1), K(m, i), ϕ, i) + R(m, i))
        end

        apply_bottom_bc!(rhs, bcs.bottom, m)
        apply_top_bc!(rhs, bcs.top, m)
    end

    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
           ϕ = solution[j]
        rhsϕ = rhs[j]
        bcsϕ = bcs[j]
          Aϕ = eqn.A[j]
          Lϕ = eqn.L[j]
          Kϕ = eqn.K[j]
          Rϕ = eqn.R[j]

        kernel!(rhsϕ, bcsϕ, ϕ, Aϕ, Lϕ, Kϕ, Rϕ, m)
    end
    return nothing
end

@inline calc_explicit_rhs!(m) = 
    calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.bcs, m.solution, m)

"Evaluate the right-hand-side of ∂ϕ∂t for an implicit time-stepping method."
function calc_implicit_rhs!(rhs, eqn, bcs, solution, m)

    N = m.grid.N

    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
           ϕ = solution[j]
        rhsϕ = rhs[j]
        bcsϕ = bcs[j]
          Aϕ = eqn.A[j]
          Kϕ = eqn.K[j]
          Rϕ = eqn.R[j]

        # Advective divergence included in lhs
        @inbounds rhsϕ[1] = ∇K∇ϕ(-zero(eltype(ϕ)), Kϕ(m, 1), ϕ, 1) + Rϕ(m, 1)

        for i = 2:N
            @inbounds rhsϕ[i] = Rϕ(m, i)
        end

        @inbounds rhsϕ[N] = (-∂zA(Aϕ(m, N+1), Aϕ(m, N), ϕ, N)
                                + ∇K∇ϕ(Kϕ(m, N+1), -zero(eltype(ϕ)), ϕ, N) + Rϕ(m, N))

        apply_bottom_bc!(rhsϕ, bcsϕ.bottom, m)
        apply_top_bc!(rhsϕ, bcsϕ.top, m)
    end
    return nothing
end

@inline calc_implicit_rhs!(m) = 
    calc_implicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.bcs, m.solution, m)

@inline K_op(m, K, iᶠ, iᶜ) = K(m, iᶠ) / Δc(m.grid, iᶠ) / Δf(m.grid, iᶜ)
@inline A_op(m, A, iᴹ, iᶜ) = A(m, iᴹ) / Δc(m.grid, iᶜ)

"Build backward Euler operator for diffusive problems."
function calc_diffusive_lhs!(lhs, eqn, solution, Δt::T, m) where T
    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
         ϕ = solution[j]
        Aϕ = eqn.A[j]
        Lϕ = eqn.L[j]
        Kϕ = eqn.K[j]
         L = lhs[j]

        for i in interiorindices(ϕ)
            @inbounds begin
                L.du[i]   = Δt * (A_op(m, Aϕ, i+1, i+1) - K_op(m, Kϕ, i+1, i))
                L.d[i]    = one(T) + Δt * (Lϕ(m, i) + K_op(m, Kϕ, i+1, i) + K_op(m, Kϕ, i, i) - A_op(m, Aϕ, i, i+1))
                L.dl[i-1] = -Δt * K_op(m, Kϕ, i, i)
            end
        end

        # Bottom row
        @inbounds L.du[1] = Δt * (A_op(m, Aϕ, 2, 2) - K_op(m, Kϕ, 2, 1))
        @inbounds L.d[1] = one(T) + Δt*(Lϕ(m, 1) + K_op(m, Kϕ, 2, 1) - A_op(m, Aϕ, 1, 2))

        # Top row
        N = length(ϕ)
        @inbounds L.dl[end] = -Δt * K_op(m, Kϕ, N, N)
        @inbounds L.d[end]  = one(T) + Δt * (Lϕ(m, N) + K_op(m, Kϕ, N, N) - A_op(m, Aϕ, N, N+1))
    end

    return nothing
end

@inline calc_diffusive_lhs!(m, Δt) =
    calc_diffusive_lhs!(m.timestepper.lhs, m.timestepper.eqn, m.solution, Δt, m)

#
# Timestepping methods
#

function update!(solution, m, eqn, bcs)
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

update!(m) = update!(m.solution, m, m.timestepper.eqn, m.bcs)

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

@inline forward_euler_step!(m, Δt) =
    forward_euler_step!(m.timestepper.rhs, m.solution, Δt)

# Forward Euler timestepping
function time_step!(m::AbstractModel{TS}, Δt) where TS <: ForwardEulerTimestepper
    update!(m)
    calc_explicit_rhs!(m)
    forward_euler_step!(m, Δt)
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

@inline backward_euler_step!(m, Δt) = 
    backward_euler_step!(m.timestepper.rhs, m.timestepper.lhs, m.solution, Δt)

"Step forward `m` by `Δt` with the backward Euler method."
function time_step!(m::AbstractModel{TS}, Δt) where TS <: BackwardEulerTimestepper
    update!(m)
    calc_implicit_rhs!(m)
    calc_diffusive_lhs!(m, Δt)
    backward_euler_step!(m, Δt)
    tick!(m.clock, Δt)
    return nothing
end
