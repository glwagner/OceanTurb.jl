module TKEMassFlux

using OceanTurb
using OceanTurb: nan2inf, inf2zero

using Printf

import ..OceanTurb: oncell, onface, maxsqrt, minuszero
import .KPP: ∂B∂z, u★, w★, isunstable
import .ModularKPP: AbstractModularKPPModel

const nsol = 5
@solution U V T S e

"Returns √ϕ if ϕ is positive and not NaN. Otherwise returns 0."
@inline function zeroed_sqrt(ϕ)
    ϕ *= 1 - isnan(ϕ)
    return maxsqrt(ϕ)
end

@inline zeroed_sqrt(ϕ, i) = @inbounds zeroed_sqrt(ϕ[i])

@inline sqrt_e(m, i) = @inbounds maxsqrt(m.solution.e[i])

@inline ∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, 
                          m.constants.g, m.constants.α, m.constants.β, i)

@inline oncell_∂B∂z(m, i) = oncell(∂B∂z, m, i) # Fallback valid for linear equations of state

@inline sqrt_∂B∂z(m, i) = maxsqrt(∂B∂z(m, i))

"Returns a velocity scale associated with convection across grid cell i."
@inline wΔ³(m, i=m.grid.N) = max(zero(eltype(m.grid)), m.state.Qb) * Δc(m.grid, i)

mutable struct Model{L, K, W, N, E, H, K0, C, ST, G, TS, S, BC, T} <: AbstractModel{TS, G, T}

                       clock :: Clock{T}
                        grid :: G
                 timestepper :: TS
                    solution :: S
                         bcs :: BC
               mixing_length :: L
          eddy_diffusivities :: K
              tke_wall_model :: W
                tke_equation :: E
        boundary_layer_depth :: H
               nonlocal_flux :: N
    background_diffusivities :: K0
                   constants :: C
                       state :: ST

end

include("state.jl")
include("nonlocal_flux.jl")
include("mixing_length.jl")
include("tke_equation.jl")
include("wall_models.jl")
include("diffusivities.jl")

function ModelBoundaryConditions(FT=Float64; U = DefaultBoundaryConditions(FT),
                                             V = DefaultBoundaryConditions(FT),
                                             T = DefaultBoundaryConditions(FT),
                                             S = DefaultBoundaryConditions(FT),
                                 )

    return (U=U, V=V, T=T, S=S)
end

function Model(; 
                                   grid = UniformGrid(N, H),
                              constants = Constants(),
                          mixing_length = SimpleMixingLength(),
                     eddy_diffusivities = IndependentDiffusivities(),
                         tke_wall_model = PrescribedSurfaceTKEFlux(),
                           tke_equation = TKEParameters(),
                   boundary_layer_depth = nothing,
                          nonlocal_flux = nothing,
               background_diffusivities = BackgroundDiffusivities(),
                                stepper = :BackwardEuler,
                                    bcs = ModelBoundaryConditions(eltype(grid))
              )

    bcs = merge(bcs, (e = TKEBoundaryConditions(eltype(grid), tke_wall_model),))

    solution = Solution((CellField(grid) for i=1:nsol)...)

    Kϕ = (U=KU, V=KV, T=KT, S=KS, e=Ke)
    Rϕ = (U=RU, V=RV, T=RT, S=RS, e=Re)
    Lϕ = (U=minuszero, V=minuszero, T=minuszero, S=minuszero, e=Le)
    eq = Equation(K=Kϕ, R=Rϕ, L=Lϕ, update=update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), 
                 grid, 
                 timestepper, 
                 solution, 
                 bcs, 
                 mixing_length, 
                 eddy_diffusivities,
                 tke_wall_model, 
                 tke_equation, 
                 boundary_layer_depth,
                 nonlocal_flux, 
                 background_diffusivities,
                 constants, 
                 State(grid, mixing_length, boundary_layer_depth, nonlocal_flux)
                )
end


@inline RU(m, i) = @inbounds   m.constants.f * m.solution.V[i] - ∂z_NLᵁ(m, i)
@inline RV(m, i) = @inbounds - m.constants.f * m.solution.U[i] - ∂z_NLⱽ(m, i)
@inline RT(m, i) = - ∂z_NLᵀ(m, i)
@inline RS(m, i) = - ∂z_NLˢ(m, i)

@inline Re(m, i) = production(m, i) + buoyancy_flux(m, i) - dissipation(m, i) - ∂z_NLᵉ(m, i)
@inline Le(m, i) = 0

end # module
