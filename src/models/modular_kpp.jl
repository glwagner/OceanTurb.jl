#= modular_kpp.jl

Here we implement a 'modular' KPP model, with three interchangable components:

1. A model for mixing depth, h
2. A model for the local diffusivity, K
3. A model for the nonlocal flux term, NL

Note below the following acronyms:

* LMD94: Large, McWilliams, and Doney (1994) "Oceanic vertical mixing..."
* RH18: Riechl and Hallberg (2018) "ePBL"
* SST07: Siebsma, Soares and Teixiera (2007) "An eddy diffusivity mass flux..."

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

using
    OceanTurb,
    LinearAlgebra

import OceanTurb: Constants

@solution U V T S

struct LMDMixingDepthParameters{T} <: AbstractParameters
    CSL  :: T  # Surface layer fraction
    CRi  :: T  # Critical bulk Richardson number
    CKE  :: T  # Unresolved turbulence parameter
end

function LMDMixingDepthParameters(T=Float64;
      CSL = 0.1,
      CRi = 0.3,
      CKE = 4.32,
     )
     LMDMixingDepthParameters{T}(CSL, CRi, CKE)
 end

struct LMDCounterGradientFluxParameters{T} <: AbstractParameters
    CNL :: T  # Non-local flux proportionality constant
end

function LMDCounterGradientFluxParameters(T=Float64; CNL=6.33)
    LMDCounterGradientFluxParameters{T}(CNL)
end

struct LMDDiffusivityParameters
    CSL   :: T  # Surface layer fraction
    Cτ    :: T  # Von Karman constant
    Cstab :: T  # Stable buoyancy flux parameter for wind-driven turbulence
    Cunst :: T  # Unstable buoyancy flux parameter for wind-driven turbulence

    Cb_U  :: T  # Buoyancy flux parameter for convective turbulence
    Cτb_U :: T  # Wind stress parameter for convective turbulence
    Cb_T  :: T  # Buoyancy flux parameter for convective turbulence
    Cτb_T :: T  # Wind stress parameter for convective turbulence

    Cd_U  :: T  # Wind mixing regime threshold for momentum
    Cd_T  :: T  # Wind mixing regime threshold for tracers

    Cn    :: T  # Exponent for effect of stable buoyancy forcing on wind mixing
    Cmτ_U :: T  # Exponent for effect of unstable buoyancy forcing on wind mixing of U
    Cmτ_T :: T  # Exponent for effect of unstable buoyancy forcing on wind mixing of T
    Cmb_U :: T  # Exponent for effect of wind on convective mixing of U
    Cmb_T :: T  # Exponent for effect of wind on convective mixing of T

    KU₀   :: T  # Interior viscosity for velocity
    KT₀   :: T  # Interior diffusivity for temperature
    KS₀   :: T  # Interior diffusivity for salinity
end

function LMDDiffusivityParameters(T=Float64;
      CSL = 0.1,
       Cτ = 0.4,
    Cstab = 2.0,
    Cunst = 6.4,
     Cb_U = 0.599,
     Cb_T = 1.36,
     Cd_U = 0.5,
     Cd_T = 2.5,
       Cn = 1.0,
    Cmτ_U = 1/4,
    Cmτ_T = 1/2,
    Cmb_U = 1/3,
    Cmb_T = 1/3,
     K₀=1e-5, KU₀=K₀, KT₀=K₀, KS₀=K₀,
     # These should not be changed under ordinary circumstances:
     Cτb_U = (Cτ / Cb_U)^(1/Cmb_U) * (1 + Cunst*Cd_U)^(Cmτ_U/Cmb_U) - Cd_U,
     Cτb_T = (Cτ / Cb_T)^(1/Cmb_T) * (1 + Cunst*Cd_T)^(Cmτ_T/Cmb_T) - Cd_T
     )

     LMDDiffusivityParameters{T}(CSL, Cτ, Cstab, Cunst,
                                 Cb_U, Cτb_U, Cb_T, Cτb_T, Cd_U, Cd_T,
                                 Cn, Cmτ_U, Cmτ_T, Cmb_U, Cmb_T,
                                 KU₀, KT₀, KS₀)
end

struct HoltslagDiffusivityParameters{T} <: AbstractParameters
    Cτ  :: T
    Cτb :: T
end

function HoltslagDiffusivityParameters(T=Float64; Cτ=0.4, Cτb=39)
    HoltslagDiffusivityParameters{T}(Cτ, Cτb)
end


struct CHCounterGradientFluxParameters{T} <: AbstractParameters
    Ca  :: T
    Cσ  :: T
    Cσw :: T
end

function CHCounterGradientFluxParameters(T=Float64;
    Ca  = 2.0,
    Cσ  = 1.3,
    Cσh = 0.6)
    CHCounterGradientFluxParameters{T}(Ca, Cσ, Cσh)
end

struct SSTMassFluxParameters{T} <: AbstractParameters
    Ce :: T
    Cμ :: T
    Cb :: T
    Cm :: T
    Cα :: T
end

function SSTMassFluxParameters(T=Float64;
    Ce = 0.4,
    Cμ = 0.15,
    Cb = 0.5,
    Cm = 0.3,
    Cα = 1.0
    )
    SSTMassFluxParameters{T}(Ce, Cμ, Cb, Cm, Cα)

# Shape functions (these shoul become parameters eventually).
# 'd' is a non-dimensional depth coordinate.
default_shape_N(d) = 0 < d < 1 ? d*(1-d)^2 : 0
const default_shape_K = default_shape_N

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

struct Model{KP, NP, HP, TS, G, T} <: AbstractModel{TS, G, T}
    @add_standard_model_fields
    K_params  :: KP
    NL_params :: NP
    h_params  :: HP
    constants :: Constants{T}
    state     :: State{T}
end

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
        h_params = LMDMixingDepthParameters(),
        K_params = LMDDiffusivityParameters(),
       NL_params = LMDCounterGradientFluxParameters(),
         stepper = :BackwardEuler,
             bcs = BoundaryConditions((ZeroFluxBoundaryConditions() for i=1:nsol)...)
    )

    solution = Solution((CellField(grid) for i=1:nsol)...)
    K = Accessory{Function}(KU, KV, KT, KS)
    R = Accessory{Function}(RU, RV, RT, RS)
    eqn = Equation(R, K, update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eqn, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, parameters, constants, State())
end


# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

## ** The K-Profile-Parameterization! **
d(m, i) = -m.grid.zf[i] / m.state.h

K_KPP(h, W, d) = 0 < d < 1 ? max(0, h * W * d * (1 - d)^2) : 0

# Holtslag K-profile
W_Holtslag(Cτ, Cτb, ωb, ωτ, h, d) = Cτ * h * ωb * ((ωτ/ωb)^3 + Cτb*d)

# Mass flux parameterzation
σw(Cσ, Cσh, ωτ, ωb, h, z) = Cσ * ωb * ((ωτ/ωb)^3 - Cσh*z/h)^(1/3) * (1 + z/h)^(1/2)
entrainment(Ce, Δz, h, z) = Ce * ( 1/(z+Δz) + 1/(h+z+Δz) )





"Return the buoyancy gradient at face point i."
∂B∂z(T, S, g, α, β, i) = g * (α*∂z(T, i) - β*∂z(S, i))
∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α, m.constants.β, i)

#
# Diagnosis of mixing depth "h" by LMD scheme
#

"Returns the surface_layer_average for mixing depth h = -zf[i]."
function surface_layer_average(c, CSL, i)
    if i > c.grid.N # Return surface value
        return onface(c, c.grid.N+1)
    else
        iε = length(c)+1 - CSL*(length(c)+1 - i) # (fractional) face "index" of the surface layer
        face = ceil(Int, iε)  # the next cell face above the fractional depth
        frac = face - iε # the fraction of the lowest cell in the surface layer.
        surface_layer_integral = convert(eltype(c), 0)

        # Contribution of fractional cell to total integral
        if frac > 0
            surface_layer_integral += frac * Δf(c, face-1) * c[face-1]
        end

        # Add cells above face, if there are any.
        for j = face:length(c)
          @inbounds surface_layer_integral += Δf(c, j) * c[j]
        end

        h = -c.grid.zf[i] # depth

        return surface_layer_integral / (CSL*h)
    end
end

"""
Return Δc(hᵢ), the difference between the surface-layer average of c and its value at depth hᵢ, where
i is a face index.
"""
Δ(c, CSL, i) = surface_layer_average(c, CSL, i) - onface(c, i)

"Returns the parameterization for unresolved KE at face point i."
function unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, β, i)
    return CKE * h^(4/3) * sqrt(max(0, Bz)) * max(0, Fb)^(1/3) + CKE₀
end

"""
    bulk_richardson_number(model, i)

Returns the bulk Richardson number of `model` at face `i`.
"""
function bulk_richardson_number(U, V, T, S, Fb, CKE, CKE₀, CSL, g, α, β, i)
    h = -U.grid.zf[i]
    # (h - hε) * ΔB
    h⁺ΔB = h * (1 - 0.5CSL) * g * (α*Δ(T, CSL, i) - β*Δ(S, CSL, i))

    Bz = ∂B∂z(T, S, g, α, β, i)
    unresolved_KE = unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, β, i)
    KE = Δ(U, CSL, i)^2 + Δ(V, CSL, i)^2 + unresolved_KE

    if KE == 0 && h⁺ΔB == 0 # Alistar Adcroft's theorem
        return 0
    else
        return h⁺ΔB / KE
    end
end

bulk_richardson_number(m, i) = bulk_richardson_number(
    m.solution.U, m.solution.V, m.solution.T, m.solution.S,
    m.state.Fb, m.parameters.CKE, m.parameters.CKE₀, m.parameters.CSL, m.constants.g,
    m.constants.α, m.constants.β, i)

"""
    mixing_depth(model)

Calculate the mixing depth 'h' for `model`.
"""
function mixing_depth(m)
    # Descend through grid until Ri rises above critical value
    ih₁ = m.grid.N + 1 # start at top.
    Ri₁ = bulk_richardson_number(m, ih₁) # should be 0.
    while ih₁ > 1 && Ri₁ < m.parameters.CRi
        ih₁ -= 1 # descend
        Ri₁ = bulk_richardson_number(m, ih₁)
    end

    # Edge cases:
    # 1. Mixing depth is 0 or whole domain:
    if ih₁ == 1 || ih₁ == m.grid.N+1
        z★ = m.grid.zf[ih₁]

    # 2. Ri is infinite somewhere inside the domain.
    elseif !isfinite(Ri₁)
        z★ = m.grid.zc[ih₁]

    # Main case: mixing depth is in the interior.
    else
        ΔRi = bulk_richardson_number(m, ih₁+1) - Ri₁ # linearly interpolate to find h.
        # x = x₀ + Δx * (y-y₀) / Δy
        z★ = m.grid.zf[ih₁] + Δf(m.grid, ih₁) * (m.parameters.CRi - Ri₁) / ΔRi
    end

    return -z★ # "depth" is negative height.
end

#
# Vertical velocity scale
#

"Return true if the boundary layer is unstable and convecting."
isunstable(model) = model.state.Fb > 0

"Return true if the boundary layer is forced."
isforced(model) = model.state.Fu != 0 || model.state.Fv != 0 || model.state.Fb != 0

"Return the turbuent velocity scale associated with wind stress."
ωτ(Fu, Fv) = (Fu^2 + Fv^2)^(1/4)
ωτ(m::AbstractModel) = ωτ(m.state.Fu, m.state.Fv)

"Return the turbuent velocity scale associated with convection."
ωb(Fb, h) = abs(h * Fb)^(1/3)
ωb(m::AbstractModel) = ωb(m.state.Fb, m.state.h)

"Return truncated, non-dimensional depth coordinate."
dϵ(m::AbstractModel, d) = min(m.parameters.CSL, d)

"Return the vertical velocity scale at depth d for a stable boundary layer."
W_KPP_stable(Cτ, Cstab, Cn, ωτ, ωb, d) = Cτ * ωτ / (1 + Cstab * d * (ωb/ωτ)^3)^Cn

"Return the vertical velocity scale at scaled depth dϵ for an unstable boundary layer."
function W_KPP_unstable(CSL, Cd, Cτ, Cunst, Cb, Cτb, Cmτ, Cmb, ωτ, ωb, d)
    dϵ = min(CSL, d)
    if dϵ < Cd * (ωτ/ωb)^3
        return Cτ * ωτ * ( 1 + Cunst * dϵ * (ωb/ωτ)^3 )^Cmτ
    else
        return Cb * ωb * ( dϵ + Cτb * (ωτ/ωb)^3 )^Cmb
    end
end

function W_KPP_unstable_U(m, i)
    return W_KPP_unstable(m.parameters.CSL, m.parameters.Cd_U, m.parameters.Cτ, m.parameters.Cunst,
                            m.parameters.Cb_U, m.parameters.Cτb_U,
                            m.parameters.Cmτ_U, m.parameters.Cmb_U,
                            ωτ(m), ωb(m), d(m, i)
                            )
end

function W_KPP_unstable_T(m, i)
    return W_KPP_unstable(m.parameters.CSL, m.parameters.Cd_T, m.parameters.Cτ, m.parameters.Cunst,
                            m.parameters.Cb_T, m.parameters.Cτb_T,
                            m.parameters.Cmτ_T, m.parameters.Cmb_T,
                            ωτ(m), ωb(m), d(m, i)
                            )
end

function W_KPP_stable(m, i)
    return W_KPP_stable(m.parameters.Cτ, m.parameters.Cstab, m.parameters.Cn,
                          ωτ(m), ωb(m), d(m, i)
                          )
end

"Return the vertical velocity scale for momentum at face point i."
function W_KPP_U(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return W_KPP_unstable_U(m, i)
    else
        return W_KPP_stable(m, i)
    end
end


"Return the vertical velocity scale for tracers at face point i."
function W_KPP_T(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return W_KPP_unstable_T(m, i)
    else
        return W_KPP_stable(m, i)
    end
end

const W_KPP_V = W_KPP_U
const W_KPP_S = W_KPP_T

#
# Non-local flux
#

"""
    N(CNL, flux, d, shape=default_shape)

Returns the nonlocal flux, N = CNL*flux*shape(d),
where `flux` is the flux of some quantity out of the surface,
`shape` is a shape function, and `d` is a non-dimensional depth coordinate
that increases from 0 at the surface to 1 at the bottom of the mixing layer.

Because flux is defined as pointing in the positive direction,
a positive surface flux implies negative surface flux divergence,
which implies a reduction to the quantity in question.
For example, positive heat flux out of the surface implies cooling.
"""
N_LMD(CNL, flux, d, shape=default_shape_N) = CNL * flux * shape(d)
N_CH(Ca, Cσ, Cσh, ωτ, ωb, Fϕ, h, z) = Ca * ωb / σw(Cσ, Cσh, ωτ, ωb, h, z)^2 * Fϕ

function ∂N∂z(CNL, Fϕ, m, i)
    if isunstable(m)
        return (N_LMD(CNL, Fϕ, d(m, i+1)) - N_LMD(CNL, Fϕ, d(m, i))) / Δf(m.grid, i)
    else
        return 0
    end
end

∂NT∂z(m, i) = ∂N∂z(m.parameters.CNL, m.state.Fθ, m, i)
∂NS∂z(m, i) = ∂N∂z(m.parameters.CNL, m.state.Fs, m, i)

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
