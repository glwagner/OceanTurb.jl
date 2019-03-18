# Timesteppers for OceanTurb.jl

"Update the clock after one iteration."
function update!(clock, Δt)
  clock.time += Δt
  clock.iter += 1
  return nothing
end

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

function unpack(model::AbstractModel{TS}, i) where TS <: ForwardEulerTimestepper
  ϕ = model.solution[i]
  rhs = model.timestepper.rhs[i]
  return ϕ, rhs
end

# Forward Euler timestepping
function iterate!(model::AbstractModel{TS}, Δt) where TS <: ForwardEulerTimestepper

  model.equation.R(model.timestepper.rhs, model)

  for j in eachindex(model.solution)
      ϕ, rhs = unpack(model, j)
      for i in eachindex(ϕ)
          @inbounds ϕ[i] += Δt * rhs[i]
      end
  end

  update!(model.clock, Δt)

  return nothing
end


#=
#
# Backward Euler Timestepper
#

struct BackwardEulerTimestepper{R, L} <: Timestepper
    rhs::R
    lhs::L
    function BackwardEulerTimestepper(solution, lhs)
        rhs = deepcopy(solution)
        for fld in rhs
            set!(fld, 0)
        end
      new{typeof(rhs), typeof(lhs)}(rhs, lhs)
    end
end

#=
function BackwardEulerTimestepper(solution)
    lhs_array = []
    for s in solution
        T = eltype(s)
        A
        N = length(s)
        lhs_s = Tridiagonal{T, A}(zeros(N-1), zeros(N), zeros(N-1))
        push!(lhs_array, lhs_s)
    end
    lhs = DiffusiveOperator(lhs_array...)

    rhs = deepcopy(solution)
    for fld in rhs
        set!(fld, 0)
    end

    BackwardEulerTimestepper(rhs, lhs)
end
=#

flux_div_op(m, κ, face, cell) = κ(m, face) / Δc(m, face) / Δf(m, cell)

# Backward Euler timestepping for problems with diffusivity
function iterate!(m::AbstractModel{TS}, Δt) where TS <: BackwardEulerTimestepper

    be = m.timestepper

    for j in eachindex(m.solution)
        ϕ, N, K, rhs, bcs = unpack(m, j)

        # Interior step
        for i in interior(ϕ)
            @inbounds begin
                # Compute rhs of backward Euler inversion equation
                rhs.data[i] = N(m, i) + ϕ.data[i]

                # Build backward euler operator for diffusive problems
                be.lhs.dl[i-1] = -Δt*flux_div_op(m, K, i, i)
                be.lhs.du[i]   = -Δt*flux_div_op(m, K, i+1, i)
                be.lhs.d[i]    = 1 + Δt*(flux_div_op(m, K, i+1, i) + flux_div_op(m, K, i, i))
            end
        end

        # Boundary conditions: inhomogeneous BCs added to rhs.
        @inbounds rhs.data[end] = N(m, bcs.top) + ϕ.data[end]
        @inbounds rhs.data[1]   = N(m, bcs.bottom) + ϕ.data[1]

        @inbounds be.lhs.du[1] = -Δt*flux_div_op(m, K, 2, 1)
        @inbounds be.lhs.d[1] = 1 + Δt*flux_div_op(m, K, 2, 1)

        @inbounds be.lhs.dl[end] = -Δt*flux_div_op(m, K, length(ϕ), length(ϕ))
        @inbounds be.lhs.d[end]  = 1 + Δt*flux_div_op(m, K, length(ϕ), length(ϕ))
    end

    # Update solution by inverting Tridiagonal lhs matrix
    for j in eachindex(m.solution)
        ϕ, N, K, rhs, bcs = unpack(m, j)
        ϕ.data .= rhs.data \ be.lhs
    end

    update!(m.clock, Δt)

  return nothing
end


function unpack(model::AbstractModel{TS}, i) where TS <: BackwardEulerTimestepper
  ϕ = model.solution[i]
  N = model.equation.N[i]
  K = model.equation.K[i]
  rhs = model.timestepper.rhs[i]
  bcs = model.bcs[i]
  return ϕ, N, K, rhs, bcs
end
=#
