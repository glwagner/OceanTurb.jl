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
    function kernel!(rhs, ϕ, M::MF, L::LF, K::KF, R::RF, m) where {MF, LF, KF, RF}
        for i in eachindex(ϕ)
            @inbounds rhs[i] = (-L(m, i) * ϕ[i] - ∂zM(M(m, i+1), M(m, i), ϕ, i)
                                    + ∇K∇ϕ(K(m, i+1), K(m, i), ϕ, i) + R(m, i))
        end
    end

    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
           ϕ = solution[j]
        rhsϕ = rhs[j]
          Mϕ = eqn.M[j]
          Lϕ = eqn.L[j]
          Kϕ = eqn.K[j]
          Rϕ = eqn.R[j]

        kernel!(rhsϕ, ϕ, Mϕ, Lϕ, Kϕ, Rϕ, m)
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
          Mϕ = eqn.M[j]
          Kϕ = eqn.K[j]
          Rϕ = eqn.R[j]

        for i in eachindex(rhsϕ)
            @inbounds rhsϕ[i] = Rϕ(m, i)
        end

        @inbounds rhsϕ[N] = (-∂zM(Mϕ(m, N+1), Mϕ(m, N), ϕ, N)
                                + ∇K∇ϕ(Kϕ(m, N+1), -zero(eltype(ϕ)), ϕ, N) + Rϕ(m, N))

        # Advective divergence included in lhs
        @inbounds rhsϕ[1] = ∇K∇ϕ(-zero(eltype(ϕ)), Kϕ(m, 1), ϕ, 1) + Rϕ(m, 1)

    end
    return nothing
end

@inline K_op(m, K, iᶠ, iᶜ) = K(m, iᶠ) / Δc(m.grid, iᶠ) / Δf(m.grid, iᶜ)
@inline M_op(m, M, iᴹ, iᶜ) = M(m, iᴹ) / Δc(m.grid, iᶜ)

"Build backward Euler operator for diffusive problems."
function calc_diffusive_lhs!(lhs, eqn, solution, cΔt::T, m) where T
    ntuple(Val(length(solution))) do j
        Base.@_inline_meta
         ϕ = solution[j]
        Mϕ = eqn.M[j]
        Lϕ = eqn.L[j]
        Kϕ = eqn.K[j]
         ℒ = lhs[j]

        for i in interiorindices(ϕ)
            @inbounds begin
                ℒ.du[i]   = cΔt * (M_op(m, Mϕ, i+1, i+1) - K_op(m, Kϕ, i+1, i))
                ℒ.d[i]    = one(T) + cΔt * (Lϕ(m, i) + K_op(m, Kϕ, i+1, i) + K_op(m, Kϕ, i, i) - M_op(m, Mϕ, i, i+1))
                ℒ.dl[i-1] = -cΔt * K_op(m, Kϕ, i, i)
            end
        end

        # Bottom row
        @inbounds ℒ.du[1] = cΔt * (M_op(m, Mϕ, 2, 2) - K_op(m, Kϕ, 2, 1))
        @inbounds ℒ.d[1] = one(T) + cΔt*(Lϕ(m, 1) + K_op(m, Kϕ, 2, 1) - M_op(m, Mϕ, 1, 2))

        # Top row
        N = length(ϕ)
        @inbounds ℒ.dl[end] = -cΔt * K_op(m, Kϕ, N, N)
        @inbounds ℒ.d[end]  = one(T) + cΔt * (Lϕ(m, N) + K_op(m, Kϕ, N, N) - M_op(m, Mϕ, N, N+1))
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
    calc_diffusive_lhs!(m.timestepper.lhs, m.timestepper.eqn, m.solution, Δt, m)
    backward_euler_step!(m.timestepper.rhs, m.timestepper.lhs, m.solution, Δt)
    tick!(m.clock, Δt)
    return nothing
end

#
# Second-order Semi-implicit Backward Difference Formula (SBDF2)
#

struct SBDF2Timestepper{E, R, L} <: Timestepper
       eqn :: E
       rhs :: R
        Rⁿ :: R
      Rⁿ⁻¹ :: R
      Φⁿ⁻¹ :: R
       lhs :: L
    function SBDF2Timestepper(eqn, solution, lhs)
        rhs = deepcopy(solution)
        Rⁿ = deepcopy(solution)
        Rⁿ⁻¹ = deepcopy(solution)
        Φⁿ⁻¹ = deepcopy(solution)
        [ set!(fld, 0) for fld in rhs   ]
        [ set!(fld, 0) for fld in Rⁿ ]
        [ set!(fld, 0) for fld in Rⁿ⁻¹ ]
        [ set!(fld, 0) for fld in Φⁿ⁻¹ ]
        new{typeof(eqn), typeof(rhs), typeof(lhs)}(eqn, rhs, Rⁿ, Rⁿ⁻¹, Φⁿ⁻¹, lhs)
    end
end

function sbdf_step!(rhs, lhs, Φⁿ, Φⁿ⁻¹, Rⁿ, Rⁿ⁻¹, Δt)
    ntuple(Val(length(Φⁿ))) do j
        Base.@_inline_meta
         lhsϕ = lhs[j]
         rhsϕ = data(  rhs[j] )
          Rϕⁿ = data(   Rⁿ[j] )
        Rϕⁿ⁻¹ = data( Rⁿ⁻¹[j] )
           ϕⁿ = data(   Φⁿ[j] )
         ϕⁿ⁻¹ = data( Φⁿ⁻¹[j] )

        for i in eachindex(ϕⁿ)
            @inbounds rhsϕ[i] = (4ϕⁿ[i] - ϕⁿ⁻¹[i] + 4Δt * Rϕⁿ[i] - 2Δt * Rϕⁿ⁻¹[i]) / 3
        end

        ldiv!(ϕⁿ, lhsϕ, rhsϕ)
    end
    return nothing
end

function store_past_data!(Rⁿ, Rⁿ⁻¹, Φⁿ, Φⁿ⁻¹)
    ntuple(Val(length(Φⁿ))) do j
        Base.@_inline_meta

          Rϕⁿ = parent(   Rⁿ[j] )
        Rϕⁿ⁻¹ = parent( Rⁿ⁻¹[j] )
           ϕⁿ = parent(   Φⁿ[j] )
         ϕⁿ⁻¹ = parent( Φⁿ⁻¹[j] )

        for i in eachindex(ϕⁿ)
            @inbounds Rϕⁿ⁻¹[i] = Rϕⁿ[i]
            @inbounds ϕⁿ⁻¹[i] = ϕⁿ[i]
        end
    end
    return nothing
end

function store_a_b!(a, b)
    ntuple(Val(length(a))) do j
        Base.@_inline_meta
         aϕ = data(a[j])
         bϕ = data(b[j])

        for i in eachindex(aϕ)
            @inbounds aϕ[i] = bϕ[i]
        end
    end
    return nothing
end

function unpack(m::AbstractModel{TS}) where TS <: SBDF2Timestepper
    return  (m.timestepper.eqn, m.solution, m.timestepper.Φⁿ⁻¹,
             m.timestepper.rhs, m.timestepper.lhs,
             m.timestepper.Rⁿ, m.timestepper.Rⁿ⁻¹)
end

"""
Step forward `m` by `Δt` using the second-order semi-implicit
backward difference formula.
"""
function iterate!(m::AbstractModel{TS}, Δt) where TS <: SBDF2Timestepper

    eqn, Φⁿ, Φⁿ⁻¹, rhs, lhs, Rⁿ, Rⁿ⁻¹ = unpack(m)

    update!(m.bcs, eqn, Φⁿ, m)
    calc_implicit_rhs!(Rⁿ, eqn, Φⁿ, m)

    #calc_diffusive_lhs!(lhs, eqn, Φⁿ, Δt, m)
    #store_a_b!(rhs, Rⁿ)
    #backward_euler_step!(rhs, lhs, Φⁿ, Δt)

    if m.clock.iter == 0 # take Backward Euler step
        calc_diffusive_lhs!(lhs, eqn, Φⁿ, Δt, m)
        store_a_b!(rhs, Rⁿ)
        backward_euler_step!(rhs, lhs, Φⁿ, Δt)
    else # sbdf step
        calc_diffusive_lhs!(lhs, eqn, Φⁿ, 2Δt/3, m)
        sbdf_step!(rhs, lhs, Φⁿ, Φⁿ⁻¹, Rⁿ, Rⁿ⁻¹, Δt)
    end

    store_past_data!(Rⁿ, Rⁿ⁻¹, Φⁿ, Φⁿ⁻¹)
    tick!(m.clock, Δt)

    return nothing
end
