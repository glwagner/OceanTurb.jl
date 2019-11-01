#= modular_kpp.jl

Here we implement a 'modular' KPP model, with three interchangable components:

1. A model for mixing depth, h
2. A model for the local diffusivity, K
3. A model for the nonlocal flux term, M

Note below the following acronyms:

* LMD94: Large, McWilliams, and Doney (1994) "Oceanic vertical mixing..."
* RH18: Riechl and Hallberg (2018) "ePBL"
* LMD07: Siebsma, Soares and Teixiera (2007) "An eddy diffusivity mass flux..."

For mixing depth models we have

1. The diagnostic bulk Richarson number scheme proposed by LMD94
2. The diagnostic TKE budget-based scheme proposed by RH18

For K-profile models we have

1. The K-profile proposed by LMD94
2. The K-profile proposed by Holtslag (1998), and described in SST07

For nonlocal flux models we have

1. The countergradient flux model proposed by LMD94
2. The mass flux model proposed by SST07

=#

module ModularKPP

export
    LMDMixingDepth,
    LMDCounterGradientFlux,
    LMDDiffusivity,
    ROMSMixingDepth,
    HoltslagDiffusivity

using
    OceanTurb,
    LinearAlgebra

import OceanTurb.KPP: 𝒲_unstable, 𝒲_stable, ωτ, ωb, d,
                      isunstable, isforced, unresolved_kinetic_energy,
                      ∂B∂z

abstract type AbstractModularKPPModel{K, H, N, TS, G, T} <: AbstractModel{TS, G, T} end

const nsol = 4
@solution U V T S

# Shape functions (these should become parameters eventually).
# 'd' is a non-dimensional depth coordinate.
default_shape_M(d) = 0 < d < 1 ? d * (1-d)^2 : 0

include("state.jl")
include("model_boundary_conditions.jl")
include("forcing.jl")

"""
    Model{S, G, T, U, B, F} <: AbstractModel{S, G, T}

Struct for KPP models.
"""
mutable struct Model{KP, NP, HP, SP, SO, BC, ST, TS, G, T, F} <: AbstractModularKPPModel{KP, NP, HP, TS, G, T}
           clock :: Clock{T}
            grid :: G
     timestepper :: TS
        solution :: SO
             bcs :: BC
     diffusivity :: KP
    nonlocalflux :: NP
     mixingdepth :: HP
        kprofile :: SP
       constants :: Constants{T}
           state :: ST
         forcing :: F
end

include("mixing_depth.jl")
include("diffusivity_profiles.jl")
include("nonlocal_flux.jl")

"""
    Model(; kwargs...)

Construct a KPP Model.
"""
function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
     diffusivity = LMDDiffusivity(),
    nonlocalflux = LMDCounterGradientFlux(),
     mixingdepth = LMDMixingDepth(),
        kprofile = DiffusivityShape(),
         stepper = :BackwardEuler,
             bcs = ModelBoundaryConditions(eltype(grid)),
         forcing = Forcing()
    )

     K = Accessory{Function}(KU, KV, KT, KS)
     R = Accessory{Function}(RU, RV, RT, RS)
    eq = Equation(K=K, R=R, update=update_state!)

       state = State(diffusivity, nonlocalflux, mixingdepth, grid)
    solution = Solution((CellField(grid) for i=1:nsol)...)
         lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs,
                 diffusivity, nonlocalflux, mixingdepth, kprofile, constants, state, 
                 forcing)
end

update_nonlocal_flux!(m) = nothing

#
# Equation specification
#

RU(m, i) =   m.constants.f * m.solution.V[i] + m.forcing.U(m, i)  
RV(m, i) = - m.constants.f * m.solution.U[i] + m.forcing.V(m, i)

# K_{U,V,T,S} is calculated at face points
KU(m::AbstractModularKPPModel{<:LMDDiffusivity}, i) =
    K_KPP(m.state.h, 𝒲_LMD_U(m, i), d(m, i), m.kprofile) + m.diffusivity.KU₀

KT(m::AbstractModularKPPModel{<:LMDDiffusivity}, i) =
    K_KPP(m.state.h, 𝒲_LMD_T(m, i), d(m, i), m.kprofile) + m.diffusivity.KT₀

KS(m::AbstractModularKPPModel{<:LMDDiffusivity}, i) =
    K_KPP(m.state.h, 𝒲_LMD_S(m, i), d(m, i), m.kprofile) + m.diffusivity.KS₀

KU(m::AbstractModularKPPModel{<:HoltslagDiffusivity}, i) =
    K_KPP(m.state.h, 𝒲_Holtslag(m, i), d(m, i), m.kprofile) + m.diffusivity.KU₀

KT(m::AbstractModularKPPModel{<:HoltslagDiffusivity}, i) =
    K_KPP(m.state.h, 𝒲_Holtslag(m, i), d(m, i), m.kprofile) + m.diffusivity.KT₀

KS(m::AbstractModularKPPModel{<:HoltslagDiffusivity}, i) =
    K_KPP(m.state.h, 𝒲_Holtslag(m, i), d(m, i), m.kprofile) + m.diffusivity.KS₀

const KV = KU

RT(m, i) = - ∂NLT∂z(m, i) + m.forcing.T(m, i)
RS(m, i) = - ∂NLS∂z(m, i) + m.forcing.S(m, i)

end # module
