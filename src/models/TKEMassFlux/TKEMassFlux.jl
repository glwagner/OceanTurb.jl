module TKEMassFlux

using OceanTurb

using Printf

import ..OceanTurb: oncell, onface
import .KPP: ∂B∂z, isunstable, ωτ, ωb
import .ModularKPP: AbstractModularKPPModel

const nsol = 5
@solution U V T S e

minuszero(args...) = -0

include("state.jl")
include("tke.jl")

struct Model{TKE, NP, H, TS, G, T, S, BC} <: AbstractModel{TS, G, T}
          clock :: Clock{T}
           grid :: G
    timestepper :: TS
       solution :: S
            bcs :: BC
            tke :: TKE
   nonlocalflux :: NP
    mixingdepth :: H
      constants :: Constants{T}
          state :: State{T}
end

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
             tke = TKEParameters(),
    nonlocalflux = nothing,
     mixingdepth = ModularKPP.LMDMixingDepth(),
         stepper = :ForwardEuler,
    )

    solution = Solution((CellField(grid) for i=1:nsol)...)

    bcs = (
        U = DefaultBoundaryConditions(eltype(grid)),
        V = DefaultBoundaryConditions(eltype(grid)),
        T = DefaultBoundaryConditions(eltype(grid)),
        S = DefaultBoundaryConditions(eltype(grid)),
        e = DefaultBoundaryConditions(eltype(grid))
    )

    Kϕ = (U=KU, V=KV, T=KT, S=KS, e=Ke)
    Rϕ = (U=RU, V=RV, T=RT, S=RS, e=Re)
    Lϕ = (U=minuszero, V=minuszero, T=minuszero, S=minuszero, e=Le)
    eq = Equation(K=Kϕ, R=Rϕ, L=Lϕ, update=update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, tke,
                 mixingdepth, constants, State())
end

@inline K(m, i) = @inbounds mixing_length(m, i) * onface(sqrt_e, m, i)

@inline KU(m, i) = m.tke.KU₀ + m.tke.CK_U * K(m, i)
@inline KT(m, i) = m.tke.KT₀ + m.tke.CK_T * K(m, i)
@inline KS(m, i) = m.tke.KS₀ + m.tke.CK_T * K(m, i)
@inline Ke(m, i) = m.tke.Ke₀ + m.tke.CK_e * K(m, i)

const KV = KU

@inline RU(m, i) = @inbounds   m.constants.f * m.solution.V[i]
@inline RV(m, i) = @inbounds - m.constants.f * m.solution.U[i]
@inline RT(m, i) = 0
@inline RS(m, i) = 0

@inline Re(m, i) = oncell(production, m, i) + oncell(buoyancy_flux, m, i) - dissipation(m, i)
@inline Le(m, i) = -0 #@inbounds m.tke.CDe / mixing_time(m, i)

end # module
