module TKEMassFlux

using OceanTurb

using Printf

import ..OceanTurb: oncell, onface
import .KPP: ∂B∂z, isunstable, ωτ, ωb
import .ModularKPP: AbstractModularKPPModel

const nsol = 5
@solution U V T S e

minuszero(args...) = -0

@inline maxsqrt(ϕ::T) where T = sqrt(max(zero(T), ϕ))
@inline maxsqrt(ϕ, i) = @inbounds sqrt(max(0, ϕ[i]))

@inline oncell(f::Function, m, i) = (f(m, i) + f(m, i+1)) / 2
@inline onface(f::Function, m, i) = (f(m, i) + f(m, i-1)) / 2

@inline sqrt_e(m, i) = @inbounds maxsqrt(m.solution.e[i])

@inline ∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α,
                          m.constants.β, i)

@inline sqrt_∂B∂z(m, i) = maxsqrt(∂B∂z(m, i))

mutable struct Model{L, P, K, TS, G, T, S, BC, C, ST} <: AbstractModel{TS, G, T}
            clock :: Clock{T}
             grid :: G
      timestepper :: TS
         solution :: S
              bcs :: BC
    mixing_length :: L
    nonlocal_flux :: P
     tke_equation :: K
        constants :: C
            state :: ST
end

include("state.jl")
include("mixing_length.jl")
include("tke_equation.jl")

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
   mixing_length = SimpleMixingLength(),
   nonlocal_flux = nothing,
    tke_equation = TKEParameters(),
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

    return Model(Clock(), grid, timestepper, solution, bcs, mixing_length, 
                 nonlocal_flux, tke_equation, constants, State())
end

@inline K(m, i) = @inbounds mixing_length_face(m, i) * onface(sqrt_e, m, i)

@inline KU(m, i) = m.tke_equation.KU₀ + m.tke_equation.CK_U * K(m, i)
@inline KT(m, i) = m.tke_equation.KT₀ + m.tke_equation.CK_T * K(m, i)
@inline KS(m, i) = m.tke_equation.KS₀ + m.tke_equation.CK_T * K(m, i)
@inline Ke(m, i) = m.tke_equation.Ke₀ + m.tke_equation.CK_e * K(m, i)

const KV = KU

@inline RU(m, i) = @inbounds   m.constants.f * m.solution.V[i]
@inline RV(m, i) = @inbounds - m.constants.f * m.solution.U[i]
@inline RT(m, i) = 0
@inline RS(m, i) = 0

@inline Re(m, i) = oncell(production, m, i) + oncell(buoyancy_flux, m, i) - dissipation(m, i)
@inline Le(m, i) = 0

end # module
