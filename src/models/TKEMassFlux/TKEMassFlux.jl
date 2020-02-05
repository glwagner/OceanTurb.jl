module TKEMassFlux

using OceanTurb
using OceanTurb: nan2inf, inf2zero

using Printf

import ..OceanTurb: oncell, onface
import .KPP: ∂B∂z, u★, isunstable
import .ModularKPP: AbstractModularKPPModel

const nsol = 5
@solution U V T S e

minuszero(args...) = -0

@inline maxsqrt(ϕ::T) where T = sqrt(max(zero(T), ϕ))
@inline maxsqrt(ϕ, i) = @inbounds sqrt(max(zero(eltype(ϕ)), ϕ[i]))

@inline oncell(f::Function, m, i) = (f(m, i) + f(m, i+1)) / 2
@inline onface(f::Function, m, i) = (f(m, i) + f(m, i-1)) / 2

@inline sqrt_e(m, i) = @inbounds maxsqrt(m.solution.e[i])

@inline ∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α,
                          m.constants.β, i)

@inline sqrt_∂B∂z(m, i) = maxsqrt(∂B∂z(m, i))

mutable struct Model{L, H, W, P, K, TS, G, T, S, BC, C, ST} <: AbstractModel{TS, G, T}
                   clock :: Clock{T}
                    grid :: G
             timestepper :: TS
                solution :: S
                     bcs :: BC
           mixing_length :: L
    boundary_layer_depth :: H
           nonlocal_flux :: P
            tke_equation :: K
          tke_wall_model :: W
               constants :: C
                   state :: ST
end

include("state.jl")
include("mixing_length.jl")
include("tke_equation.jl")

function Model(; 
                      grid = UniformGrid(N, L),
                 constants = Constants(),
             mixing_length = SimpleMixingLength(),
      boundary_layer_depth = nothing,
             nonlocal_flux = nothing,
              tke_equation = TKEParameters(),
            tke_wall_model = nothing,
                   stepper = :ForwardEuler,
)

    solution = Solution((CellField(grid) for i=1:nsol)...)

    tke_bcs = TurbulentKineticEnergyBoundaryConditions(eltype(grid), tke_wall_model)

    bcs = (
        U = DefaultBoundaryConditions(eltype(grid)),
        V = DefaultBoundaryConditions(eltype(grid)),
        T = DefaultBoundaryConditions(eltype(grid)),
        S = DefaultBoundaryConditions(eltype(grid)),
        e = tke_bcs
    )

    Kϕ = (U=KU, V=KV, T=KT, S=KS, e=Ke)
    Rϕ = (U=RU, V=RV, T=RT, S=RS, e=Re)
    Lϕ = (U=minuszero, V=minuszero, T=minuszero, S=minuszero, e=Le)
    eq = Equation(K=Kϕ, R=Rϕ, L=Lϕ, update=update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, mixing_length, boundary_layer_depth,
                 nonlocal_flux, tke_equation, tke_wall_model, constants, 
                 State(mixing_length, boundary_layer_depth))
end

@inline K(m, i) = @inbounds diffusivity_mixing_length(m, i) * onface(sqrt_e, m, i)

@inline KU(m, i) = m.tke_equation.KU₀ + m.tke_equation.Cᴷᵤ * K(m, i)
@inline KT(m, i) = m.tke_equation.KT₀ + m.tke_equation.Cᴷᵩ * K(m, i)
@inline KS(m, i) = m.tke_equation.KS₀ + m.tke_equation.Cᴷᵩ * K(m, i)
@inline Ke(m, i) = m.tke_equation.Ke₀ + m.tke_equation.Cᴷₑ * K(m, i)

const KV = KU

@inline RU(m, i) = @inbounds   m.constants.f * m.solution.V[i]
@inline RV(m, i) = @inbounds - m.constants.f * m.solution.U[i]
@inline RT(m, i) = 0
@inline RS(m, i) = 0

@inline Re(m, i) = oncell(production, m, i) + oncell(buoyancy_flux, m, i) - dissipation(m, i)
@inline Le(m, i) = 0

end # module
