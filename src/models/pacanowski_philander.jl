module PacanowskiPhilander

using
    StaticArrays,
    OceanTurb

export
    Parameters,
    Model

import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top, Constants

import .KPP: ∂B∂z

const nsol = 4
@solution U V T S

struct Parameters{T} <: AbstractParameters
    Cν₀ :: T
    Cν₁ :: T
    Cκ₀ :: T
    Cκ₁ :: T
    Cc  :: T
    Cn  :: T
end

function Parameters(T=Float64;
    ν₀ = 1e-4,
    ν₁ = 1e-2,
    κ₀ = 1e-5,
    κ₁ = 1e-2,
    c  = 5,
    n  = 2
    )

    Parameters{T}(ν₀, ν₁, κ₀, κ₁, c, n)
end

struct Model{TS, G, T} <: AbstractModel{TS, G, T}
    @add_standard_model_fields
    parameters :: Parameters{T}
    constants  :: Constants{T}
end

function Model(; N=10, L=1.0,
          grid = UniformGrid(N, L),
     constants = Constants(),
    parameters = Parameters(),
       stepper = :ForwardEuler,
           bcs = BoundaryConditions((ZeroFluxBoundaryConditions() for i=1:nsol)...)
    )

    solution = Solution((CellField(grid) for i=1:nsol)...)
    K = Accessory{Function}(KU, KV, KT, KS)
    R = Accessory{Any}(RU, RV, nothing, nothing)
    eqn = Equation(R, K)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eqn, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, parameters, constants)
end

#
# Equation specification
#

function local_richardson(U, V, T, S, g, α, β, i)
    Bz = ∂B∂z(T, S, g, α, β, i)
    S² = ∂z(U, i)^2 + ∂z(V, i)^2

    if S² == 0 && Bz == 0 # Alistair Adcroft's theorem
        return 0
    else
        return Bz / S²
    end
end

local_richardson(m, i) = local_richardson(m.solution.U, m.solution.V, m.solution.T,
                                          m.solution.S, m.constants.g, m.constants.α,
                                          m.constants.β, i)

#
# Equation specification
#

KU(Ri, ν₀, ν₁, c, n) = ν₀ + ν₁ / (1 + c*Ri)^n
KT(Ri, κ₀, κ₁, c, n) = κ₀ + κ₁ / (1 + c*Ri)^(n+1)

KU(m, i) = KU(local_richardson(m, i), m.parameters.Cν₀, m.parameters.Cν₁,
              m.parameters.Cc, m.parameters.Cn)

KT(m, i) = KT(local_richardson(m, i), m.parameters.Cκ₀, m.parameters.Cκ₁,
              m.parameters.Cc, m.parameters.Cn)

const KV = KU
const KS = KT

RU(m, i) =   m.constants.f * m.solution.V[i]
RV(m, i) = - m.constants.f * m.solution.U[i]

end # module
