module Diffusion

using
    StaticArrays,
    OceanTurb

export
    Parameters,
    Model

import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top

import .KPP: Constants, ∂B∂z

const nsol = 4

# Just one field: "c"
@specify_solution CellField U V T S

struct Parameters{T} <: AbstractParameters
    Cν₀ :: T
    Cν₁ :: T
    Cκ₀ :: T
    Cκ₁ :: T
    Cc  :: T
    Cn  :: T
end

function Parameters(T=Float64;
    Cν₀ = 1e-4,
    Cν₁ = 1e-2,
    Cκ₀ = 1e-5,
    Cκ₁ = 1e-2,
    Cc  = 5
    Cn  = 2
    )

    Parameters{T}(Cν₀, Cν₁, Cκ₀, Cκ₁, Cc, Cn)
end

mutable struct State{T} <: FieldVector{4, T}
    Fu :: T
    Fv :: T
    Fθ :: T
    Fs :: T
end

State(T=Float64) = State{T}(0, 0, 0, 0)

"""
    update_state!(model)

Update the top flux conditions and mixing depth for `model`
and store in `model.state`.
"""
function update_state!(m)
    m.state.Fu = getbc(m, m.bcs.U.top)
    m.state.Fv = getbc(m, m.bcs.V.top)
    m.state.Fθ = getbc(m, m.bcs.T.top)
    m.state.Fs = getbc(m, m.bcs.S.top)
    return nothing
end

struct Model{P, TS, G, E, T} <: AbstractModel{TS, G, E, T}
    @add_standard_model_fields
    parameters :: Parameters{T}
    constants  :: Constants{T}
    state      :: State{T}
end

function Model(; N=10, L=1.0,
          grid = UniformGrid(N, L),
     constants = Constants(),
    parameters = Parameters(),
       stepper = :ForwardEuler,
           bcs = BoundaryConditions((ZeroFluxBoundaryConditions() for i=1:nsol)...)
    )

    solution = Solution((CellField(grid) for i=1:nsol)...)
    equation = Equation(calc_rhs_explicit!)
    timestepper = Timestepper(:ForwardEuler, solution)

    return Model(timestepper, grid, equation, solution, bcs, Clock(),
                    parameters, constants, State())
end

#
# Equation specification
#

function local_richardson(U, V, T, S, g, α, β, i)
    return ∂B∂z(T, S, g, α, β, i) / (∂z(U, i)^2 + ∂z(V, i)^2)
end

local_richardson(m, i) = local_richardson(m.solution.U, m.solution.V, m.solution.T,
                                          m.solution.S, m.constants.g, m.constants.α, m.constants.β)

# Equation specification
KU(Ri, ν₀, ν₁, c, n) = ν₀ + ν₁ / (1 + c*Ri)^n
KT(Ri, ν₀, ν₁, c, n) = κ₀ + κ₁ / (1 + c*Ri)^(n+1)

KU(m, i) = KU(local_richardon(m, i), m.parameters.Cν₀, m.parameters.Cν₁,
              m.parameters.Cc, m.parameters.Cn)

KT(m, i) = KT(local_richardon(m, i), m.parameters.Cκ₀, m.parameters.Cκ₁,
              m.parameters.Cc, m.parameters.Cn)

const KV = KU
const KS = KT

function calc_rhs_explicit!(∂t, m)

    update_state!(m)
    U, V, T, S = m.solution

    @inbounds begin

        for i in interior(∂t)
            ∂t.U[i] = ∇K∇c(KU(m, i+1), KU(m, i), U, i) + m.constants.f * V[i]
            ∂t.V[i] = ∇K∇c(KV(m, i+1), KV(m, i), V, i) - m.constants.f * U[i]
            ∂t.T[i] = ∇K∇c(KT(m, i+1), KT(m, i), T, i)
            ∂t.S[i] = ∇K∇c(KS(m, i+1), KS(m, i), S, i)
        end

        # Flux into the top (the only boundary condition allowed)
        i = m.grid.N
        ∂t.U[i] = ∇K∇c_top(KU(m, i), U, m.state.Fu) + m.constants.f*V[i]
        ∂t.V[i] = ∇K∇c_top(KV(m, i), V, m.state.Fv) - m.constants.f*U[i]
        ∂t.T[i] = ∇K∇c_top(KT(m, i), T, m.state.Fθ)
        ∂t.S[i] = ∇K∇c_top(KS(m, i), S, m.state.Fs)

        # Bottom
        i = 1
        ∂t.U[i] = ∇K∇c_bottom(KU(m, i+1), KU(m, i), U, m.bcs.U.bottom, m) + m.constants.f*V[i]
        ∂t.V[i] = ∇K∇c_bottom(KV(m, i+1), KV(m, i), V, m.bcs.V.bottom, m) - m.constants.f*U[i]
        ∂t.T[i] = ∇K∇c_bottom(KT(m, i+1), KT(m, i), T, m.bcs.T.bottom, m)
        ∂t.S[i] = ∇K∇c_bottom(KS(m, i+1), KS(m, i), S, m.bcs.S.bottom, m)
    end

    return nothing
end

end # module
