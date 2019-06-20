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
2. The K-profile proposed by Holtslag ??

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

Base.@kwdef struct LMDMixingDepth{T} <: AbstractParameters
     CSL :: T = 0.1   # Surface layer fraction
     CRi :: T = 0.3   # Critical bulk Richardson number
     CKE :: T = 4.32  # Unresolved turbulence parameter
    CKE₀ :: T = 1e-11 # Minimum unresolved turbulence kinetic energy
end

Base.@kwdef struct ROMSMixingDepth{T} <: AbstractParameters
     CSL :: T = 0.1  # Surface layer fraction
     CRi :: T = 0.3  # Critical bulk Richardson number
     CKE :: T = 5.07 # Minimum unresolved turbulence kinetic energy
     CEk :: T = 0.0  # Turbulent Ekman depth parameter
end

Base.@kwdef struct LMDCounterGradientFlux{T} <: AbstractParameters
    CNL :: T = 6.33 # Mass flux proportionality constant
end

Base.@kwdef struct DiffusivityShape{T} <: AbstractParameters
    CS0 :: T = 0.0
    CS1 :: T = 1.0
end

Base.@kwdef struct LMDDiffusivity{T} <: AbstractParameters
     CKSL :: T = 0.1   # Surface layer fraction
       Cτ :: T = 0.4   # Von Karman constant

    Cstab :: T = 2.0   # Stable buoyancy flux parameter for wind-driven turbulence
    Cunst :: T = 6.4   # Unstable buoyancy flux parameter for wind-driven turbulence

       Cn :: T = 1.0   # Exponent for effect of stable buoyancy forcing on wind mixing
    Cmτ_U :: T = 0.25  # Exponent for effect of unstable buoyancy forcing on wind mixing of U
    Cmτ_T :: T = 0.5   # Exponent for effect of unstable buoyancy forcing on wind mixing of T
    Cmb_U :: T = 1/3   # Exponent for the effect of wind on convective mixing of U
    Cmb_T :: T = 1/3   # Exponent for effect of wind on convective mixing of T

     Cd_U :: T = 0.5   # Wind mixing regime threshold for momentum
     Cd_T :: T = 2.5   # Wind mixing regime threshold for tracers

     Cb_U :: T = 0.599 # Buoyancy flux parameter for convective turbulence
     Cb_T :: T = 1.36  # Buoyancy flux parameter for convective turbulence
    Cτb_U :: T = (Cτ / Cb_U)^(1/Cmb_U) * (1 + Cunst*Cd_U)^(Cmτ_U/Cmb_U) - Cd_U  # Wind stress parameter for convective turbulence
    Cτb_T :: T = (Cτ / Cb_T)^(1/Cmb_T) * (1 + Cunst*Cd_T)^(Cmτ_T/Cmb_T) - Cd_T  # Wind stress parameter for convective turbulence

      KU₀ :: T = 1e-6 # Interior viscosity for velocity
      KT₀ :: T = 1e-7 # Interior diffusivity for temperature
      KS₀ :: T = 1e-9 # Interior diffusivity for salinity
end

Base.@kwdef struct HoltslagDiffusivity{T} <: AbstractParameters
     Cτ :: T = 0.4
    Cτb :: T = 15.6
    KU₀ :: T = 1e-6 # Interior viscosity for velocity
    KT₀ :: T = 1e-7 # Interior diffusivity for temperature
    KS₀ :: T = 1e-9 # Interior diffusivity for salinity
end

Base.@kwdef struct BulkPlumeParameters{T} <: AbstractParameters
     Ce :: T = 0.4
     Cμ :: T = 0.15
     Cb :: T = 0.5
     Cm :: T = 0.3
     Cα :: T = 1.0
     Cσ :: T = 1.0
    Cσb :: T = 1.0
end

mutable struct State{T, H, U, W}
          Fu :: T
          Fv :: T
          Fθ :: T
          Fs :: T
          Fb :: T
           h :: T
      h_crit :: H
     plume_T :: U
     plume_S :: U
    plume_w² :: W
end

plumes(args...) = nothing, nothing, nothing
h_criterion(args...) = nothing
h_criterion(::ROMSMixingDepth, grid) = FaceField(grid)

function State(diffusivity, nonlocalflux, mixingdepth, grid, T=Float64)
    plume_T, plume_S, plume_w² = plumes(nonlocalflux, grid)
    h_crit = h_criterion(mixingdepth, grid)
    State(zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
            h_crit, plume_T, plume_S, plume_w²)
end

mutable struct Model{KP, NP, HP, SP, SO, BC, ST, TS, G, T} <: AbstractModularKPPModel{KP, NP, HP, TS, G, T}
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
end

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
     diffusivity = LMDDiffusivity(),
    nonlocalflux = LMDCounterGradientFlux(),
     mixingdepth = LMDMixingDepth(),
        kprofile = DiffusivityShape(),
         stepper = :BackwardEuler
    )

     K = Accessory{Function}(KU, KV, KT, KS)
     R = Accessory{Function}(RU, RV, RT, RS)
    eq = Equation(K=K, R=R, update=update_state!)

    bcs = (
        U = DefaultBoundaryConditions(eltype(grid)),
        V = DefaultBoundaryConditions(eltype(grid)),
        T = DefaultBoundaryConditions(eltype(grid)),
        S = DefaultBoundaryConditions(eltype(grid))
    )

       state = State(diffusivity, nonlocalflux, mixingdepth, grid)
    solution = Solution((CellField(grid) for i=1:nsol)...)
         lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs,
                 diffusivity, nonlocalflux, mixingdepth, kprofile, constants, state)
end

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
    update_mixing_depth!(m)
    update_nonlocal_flux!(m)
    return nothing
end

function update_mixing_depth!(m::Model{K, NL, <:LMDMixingDepth}) where {K, NL}
    m.state.h  = mixing_depth(m)
    return nothing
end

h_weight(h, CSL, zf, i) = @inbounds -zf[i] / (CSL*h - zf[i])
h_weight(m, i) = h_weight(m.state.h, m.mixingdepth.CSL, m.grid.zf, i)

function h_kernel(U, V, T, S, CRi, CEk, g, α, β, f, i)
    @inbounds ∂z(U, i)^2 + ∂z(V, i)^2 - ∂B∂z(T, S, g, α, β, i)/CRi - CEk*f^2
end

h_kernel(m, i) = h_kernel(m.solution.U, m.solution.V, m.solution.T, m.solution.S,
                            m.mixingdepth.CRi, m.mixingdepth.CEk,
                            m.constants.g, m.constants.α, m.constants.β, m.constants.f, i)

function unresolved_kinetic_energy(m, i)
    @inbounds unresolved_kinetic_energy(-m.grid.zf[i],
        ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α, m.constants.β, i),
        m.state.Fb, m.mixingdepth.CKE, 0, m.constants.g, m.constants.α, m.constants.β)
end

"Calculate the mixing depth criterion function by integrating from z=0 downwards."
function mixing_depth_criterion!(h_crit, m)
    @inbounds h_crit[m.grid.N+1] = 0

    for i = m.grid.N:-1:1
        @inbounds h_crit[i] = h_crit[i+1] + h_weight(m, i) * h_kernel(m, i) * Δc(m.grid, i)
    end

    for i in eachindex(h_crit)
        @inbounds h_crit[i] -= unresolved_kinetic_energy(m, i) / m.grid.zf[i]
    end

    return nothing
end

linear_interp(y★, x₀, y₀, Δx, Δy) = x₀ + Δx * (y★ - y₀) / Δy

function mixing_depth(m::Model{K, NL, <:ROMSMixingDepth}) where {K, NL}
    ih₁ = findprev(x -> x<=0, m.state.h_crit.data, m.grid.N)
    @inbounds begin
        if ih₁ === nothing # Mixing depth is entire grid
            z★ = m.grid.zf[1]
        elseif ih₁ == m.grid.N # Mixing depth at surface?
            z★ = ifelse(m.state.h_crit[ih₁]==0, m.grid.zf[m.grid.N], m.grid.zf[m.grid.N+1])
        else # linearly interpolate
            # x = x₀ + Δx * (y-y₀) / Δy
            z★ = linear_interp(0, m.grid.zf[ih₁], m.state.h_crit[ih₁], Δf(m.grid, ih₁),
                                m.state.h_crit[ih₁+1] - m.state.h_crit[ih₁])
        end
    end

    return -z★
end

function update_mixing_depth!(m::Model{K, NL, <:ROMSMixingDepth}) where {K, NL}
    mixing_depth_criterion!(m.state.h_crit, m)
    m.state.h = mixing_depth(m)
    return nothing
end

update_nonlocal_flux!(m) = nothing


#
# Mixing depth
#

bulk_richardson_number(m::AbstractModel, i) = KPP.bulk_richardson_number(
    m.solution.U, m.solution.V, m.solution.T, m.solution.S,
    m.state.Fb, m.mixingdepth.CKE, m.mixingdepth.CKE₀, m.mixingdepth.CSL, m.constants.g,
    m.constants.α, m.constants.β, i)

"""
    mixing_depth(model)

Calculate the mixing depth 'h' for `model`.
"""
function mixing_depth(m)
    ih₁ = m.grid.N + 1 # start at top.
    Ri₁ = bulk_richardson_number(m, ih₁) # should be 0.

    # Descend through grid until Ri rises above critical value
    while ih₁ > 1 && Ri₁ < m.mixingdepth.CRi
        ih₁ -= 1 # descend
        Ri₁ = bulk_richardson_number(m, ih₁)
    end

    # Edge cases:
    # 1. Mixing depth is 0:
    if ih₁ == m.grid.N + 1
        z★ = m.grid.zf[ih₁]

    # 2. Mixing depth is whole domain because Ri is always less than CRi:
    elseif ih₁ == 1 && Ri₁ < m.mixingdepth.CRi
        z★ = m.grid.zf[ih₁]

    # 3. Ri is infinite somewhere inside the domain.
    elseif !isfinite(Ri₁)
        z★ = m.grid.zc[ih₁]

    # Main case: mixing depth is in the interior.
    else # Ri₁ > CRi
        ΔRi = bulk_richardson_number(m, ih₁+1) - Ri₁ # <0 linearly interpolate to find h.
        # x = x₀ + Δx * (y-y₀) / Δy
        z★ = m.grid.zf[ih₁] + Δf(m.grid, ih₁) * (m.mixingdepth.CRi - Ri₁) / ΔRi
    end

    -z★ < 0 && @warn "mixing depth $(-z★) is negative"

    return -z★ # "depth" is negative height.
end

#
# Diffusivity
#

k_profile(d, p::DiffusivityShape) = d * (1-d) * ( p.CS0 + p.CS1*(1-d) )

## ** The K-Profile-Parameterization **
K_KPP(h, 𝒲, d::T, p) where T = 0<d<1 ? max(zero(T), h * 𝒲 * k_profile(d, p)) : -zero(T)

𝒲_Holtslag(Cτ, Cτb, ωτ, ωb, d) = Cτ * (ωτ^3 + Cτb * d * ωb^3)^(1/3)
𝒲_Holtslag(m, i) = 𝒲_Holtslag(m.diffusivity.Cτ, m.diffusivity.Cτb, KPP.ωτ(m), KPP.ωb(m), KPP.d(m, i))

𝒲_LMD_unstable_U(m, i) = KPP.𝒲_unstable(
    m.diffusivity.CKSL, m.diffusivity.Cd_U,
    m.diffusivity.Cτ, m.diffusivity.Cunst,
    m.diffusivity.Cb_U, m.diffusivity.Cτb_U,
    m.diffusivity.Cmτ_U, m.diffusivity.Cmb_U,
    ωτ(m), ωb(m), d(m, i)
    )

𝒲_LMD_unstable_T(m, i) = KPP.𝒲_unstable(
    m.diffusivity.CKSL, m.diffusivity.Cd_T,
    m.diffusivity.Cτ, m.diffusivity.Cunst,
    m.diffusivity.Cb_T, m.diffusivity.Cτb_T,
    m.diffusivity.Cmτ_T, m.diffusivity.Cmb_T,
    ωτ(m), ωb(m), d(m, i)
    )

𝒲_LMD_stable(m, i) = KPP.𝒲_stable(
    m.diffusivity.Cτ, m.diffusivity.Cstab, m.diffusivity.Cn,
    ωτ(m), ωb(m), d(m, i)
    )

"Return the vertical velocity scale for momentum at face point i"
function 𝒲_LMD_U(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return 𝒲_LMD_unstable_U(m, i)
    else
        return 𝒲_LMD_stable(m, i)
    end
end

"Return the vertical velocity scale for tracers at face point i."
function 𝒲_LMD_T(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return 𝒲_LMD_unstable_T(m, i)
    else
        return 𝒲_LMD_stable(m, i)
    end
end

const 𝒲_LMD_V = 𝒲_LMD_U
const 𝒲_LMD_S = 𝒲_LMD_T

#
# Mass flux
#

# Shape functions (these shoul become parameters eventually).
# 'd' is a non-dimensional depth coordinate.
default_shape_M(d) = 0 < d < 1 ? d * (1-d)^2 : 0


function ∂NLT∂z(m::Model{K, <:LMDCounterGradientFlux}, i) where K
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Fθ, d(m, i+1), d(m, i), Δf(m.grid, i), m)
end

function ∂NLS∂z(m::Model{K, <:LMDCounterGradientFlux}, i) where K
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Fs, d(m, i+1), d(m, i), Δf(m.grid, i), m)
end

σw(ωb, ωτ, Cσ, Cσb, d) = Cσ * (ωτ^3 + Cσb * ωb^3 * d)^(1/3) * (1 - d)^(1/2)

entrainment(Ce, h, Δz, z) = Ce * (- 1 / (z + Δz) + 1 / (h + z + Δz))

function plume_buoyancy(plume_T, plume_S, T, S, α, β, g, i)
    @inbounds g*(α*(plume_T[i] - T[i]) - β*(plume_S[i] - S[i]))
end

#
# Equation specification
#

RU(m, i) =   m.constants.f * m.solution.V[i]
RV(m, i) = - m.constants.f * m.solution.U[i]

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

RT(m, i) = - ∂NLT∂z(m, i)
RS(m, i) = - ∂NLS∂z(m, i)

end # module
