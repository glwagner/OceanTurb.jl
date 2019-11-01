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
    HoltslagDiffusivity,
    DiagnosticPlumeModel

using
    OceanTurb,
    LinearAlgebra

using OceanTurb: ∂z⁺ # upwards-biased difference

import OceanTurb.KPP: 𝒲_unstable, 𝒲_stable, ωτ, ωb, d,
                      isunstable, isforced, unresolved_kinetic_energy,
                      ∂B∂z

abstract type AbstractModularKPPModel{K, H, N, TS, G, T} <: AbstractModel{TS, G, T} end

const nsol = 4
@solution U V T S

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

include("state.jl")
include("model_boundary_conditions.jl")
include("forcing.jl")
include("mixing_depth.jl")
include("diffusivity_models.jl")
include("shape_functions.jl")
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
        kprofile = StandardCubicPolynomial(),
         stepper = :BackwardEuler,
             bcs = ModelBoundaryConditions(eltype(grid)),
         forcing = Forcing()
    )

     K = Accessory{Function}(KU, KV, KT, KS)
     R = Accessory{Function}(RU, RV, RT, RS)
     M = Accessory{Function}(MU, MV, MT, MS)
    eq = Equation(K=K, R=R, M=M, update=update_state!)

       state = State(diffusivity, nonlocalflux, mixingdepth, grid)
    solution = Solution((CellField(grid) for i=1:nsol)...)
         lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs,
                 diffusivity, nonlocalflux, mixingdepth, kprofile, constants, state, 
                 forcing)
end

#
# Equation specification
#

RU(m, i) =   m.constants.f * m.solution.V[i] + m.forcing.U(m, i)  
RV(m, i) = - m.constants.f * m.solution.U[i] + m.forcing.V(m, i)

const LMDModel = AbstractModularKPPModel{<:LMDDiffusivity}
const HoltslagModel = AbstractModularKPPModel{<:HoltslagDiffusivity}

# K_{U,V,T,S} is calculated at face points
KU(m::LMDModel, i) = K_KPP(m.state.h, 𝒲_LMD_U(m, i), d(m, i), m.kprofile) + m.diffusivity.KU₀
KT(m::LMDModel, i) = K_KPP(m.state.h, 𝒲_LMD_T(m, i), d(m, i), m.kprofile) + m.diffusivity.KT₀
KS(m::LMDModel, i) = K_KPP(m.state.h, 𝒲_LMD_S(m, i), d(m, i), m.kprofile) + m.diffusivity.KS₀

KU(m::HoltslagModel, i) = K_KPP(m.state.h, 𝒲_Holtslag(m, i), d(m, i), m.kprofile) + m.diffusivity.KU₀
KT(m::HoltslagModel, i) = K_KPP(m.state.h, 𝒲_Holtslag(m, i), d(m, i), m.kprofile) + m.diffusivity.KT₀
KS(m::HoltslagModel, i) = K_KPP(m.state.h, 𝒲_Holtslag(m, i), d(m, i), m.kprofile) + m.diffusivity.KS₀

const KV = KU

RT(m, i) = - ∂z_explicit_nonlocal_flux_T(m, i) + m.forcing.T(m, i)
RS(m, i) = - ∂z_explicit_nonlocal_flux_S(m, i) + m.forcing.S(m, i)

MU(m, i) = 0
MV(m, i) = 0
MT(m, i) = mass_flux(m, i)
MS(m, i) = mass_flux(m, i)

end # module
