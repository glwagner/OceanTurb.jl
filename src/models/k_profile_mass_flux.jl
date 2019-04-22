module KProfileMassFlux

using
    OceanTurb,
    LinearAlgebra

import OceanTurb: Constants

@solution U V T S

struct KProfileParameters{T} <: AbstractParameters
    Cτ  :: T
    Cτb :: T
end

KProfileParameters(T=Float64; Cτ=0.4, Cτb=39) = KProfileParameters{T}(Cτ, Cτb)

struct MassFluxParameters{T} <: AbstractParameters
    Ce :: T
    Cμ :: T
    Cb :: T
    Cm :: T
    Cα :: T
end

function MassFluxParameters(T=Float64;
    Ce = 0.4,
    Cμ = 0.15,
    Cb = 0.5,
    Cm = 0.3,
    Cα = 1.0
    )
    MassFluxParameters{T}(Ce, Cμ, Cb, Cm, Cα)
end

mutable struct State{T} <: FieldVector{6, T}
    Fu :: T
    Fv :: T
    Fθ :: T
    Fs :: T
    Fb :: T
    h  :: T
end

State(T=Float64) = State{T}(0, 0, 0, 0, 0, 0)

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
    m.state.Fb = m.constants.g * (m.constants.α * m.state.Fθ - m.constants.β * m.state.Fs)
    m.state.h  = mixing_depth(m)
    return nothing
end

struct Model{TS, G, T} <: AbstractModel{TS, G, T}
    @add_standard_model_fields
    kparams  :: KProfileParameters{T}
    mfparams :: MassFluxParameters{T}
    constants :: Constants{T}
    state     :: State{T}
end

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
         kparams = KProfileParameters(),
        mfparams = MassFluxParameters(),
         stepper = :BackwardEuler,
             bcs = BoundaryConditions((ZeroFluxBoundaryConditions() for i=1:nsol)...)
    )

    solution = Solution((CellField(grid) for i=1:nsol)...)
    K = Accessory{Function}(KU, KV, KT, KS)
    R = Accessory{Function}(RU, RV, RT, RS)
    eqn = Equation(R, K, update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eqn, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, kparams, mfparams, State())
end

#
# Equation specification
#
# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

## ** The K-Profile-Parameterization! **
d(m, i) = -m.grid.zf[i] / m.state.h

K_KPP(h, W, d) = 0 < d < 1 ? max(0, h * W * d * (1 - d)^2) : 0

# Holtslag K-profile
W_Holtslag(Cτ, Cτb, ωb, ωτ, h, d) = Cτ * h * (ωτ^3 + Cτb * ωb^3 * d)
W_Holtslag(m, i) = W_Holtslag(m.kparams.Cτ, m.kparams.Cτb, ωb(m), ωτ(m), model.state.h, d(m, i))

# Mass flux parameterzation
σw(Cσ, Cσh, ωτ, ωb, d) = Cσ * (ωτ^3 - Cσh*ωb^3*d)^(1/3) * sqrt(1 + d)
entrainment(Ce, Δz, h, z) = Ce * (1/(z + Δz) + 1/(h + z + Δz))

"Return the buoyancy gradient at face point i."
∂B∂z(T, S, g, α, β, i) = g * (α*∂z(T, i) - β*∂z(S, i))
∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α, m.constants.β, i)

#
# Equation specification
#

RU(m, i) =   m.constants.f * m.solution.V[i]
RV(m, i) = - m.constants.f * m.solution.U[i]

# K_{U,V,T,S} is calculated at face points
KU(m::Model{<:LMDDiffusivityParameters}, i) = K_KPP(m.state.h, W_KPP_U(m, i), d(m, i)) + m.parameters.KU₀
KT(m::Model{<:LMDDiffusivityParameters}, i) = K_KPP(m.state.h, W_KPP_T(m, i), d(m, i)) + m.parameters.KT₀
KS(m::Model{<:LMDDiffusivityParameters}, i) = K_KPP(m.state.h, W_KPP_S(m, i), d(m, i)) + m.parameters.KS₀

KU(m::Model{<:HoltslagDiffusivityParameters}, i) = K_KPP(m.state.h, W_Holtslag(m, i), d(m, i)) + m.parameters.KU₀
KT(m::Model{<:HoltslagDiffusivityParameters}, i) = K_KPP(m.state.h, W_Holtslag(m, i), d(m, i)) + m.parameters.KT₀
KS(m::Model{<:HoltslagDiffusivityParameters}, i) = K_KPP(m.state.h, W_Holtslag(m, i), d(m, i)) + m.parameters.KS₀

const KV = KU

RT(m, i) = - ∂NT∂z(m, i)
RS(m, i) = - ∂NS∂z(m, i)

end # module
