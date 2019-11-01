module PacanowskiPhilander

using OceanTurb

import .KPP: ∂B∂z

const nsol = 4
@solution U V T S

Base.@kwdef struct Parameters{T} <: AbstractParameters
    Cν₀ :: T = 1e-4
    Cν₁ :: T = 1e-2
    Cκ₀ :: T = 1e-5
    Cκ₁ :: T = 1e-2
    Cc  :: T = 5.0
    Cn  :: T = 2.0
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
    K = (U=KU, V=KV, T=KT, S=KS)
    R = (U=RU, V=RV, T=RT, S=RS)
    eqn = Equation(K=K, R=R)
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
RT(m, i) = 0
RS(m, i) = 0

end # module
