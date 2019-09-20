module KPP

using OceanTurb

import OceanTurb: Constants

using Base: @propagate_inbounds

const nsol = 4
@solution U V T S

"""
    Parameters(; kwargs...)

Construct KPP parameters.
"""
@Base.kwdef struct Parameters{T<:AbstractFloat} <: AbstractParameters
    CSL   :: T  = 0.1   # Surface layer fraction
    Cτ    :: T  = 0.4   # Von Karman constant
    CNL   :: T  = 6.33  # Non-local flux proportionality constant

    Cstab :: T  = 2.0   # Stable buoyancy flux parameter for wind-driven turbulence
    Cunst :: T  = 6.4   # Unstable buoyancy flux parameter for wind-driven turbulence

       Cn :: T  = 1.0   # Exponent for effect of stable buoyancy forcing on wind mixing
    Cmτ_U :: T  = 0.25  # Exponent for effect of unstable buoyancy forcing on wind mixing of U
    Cmτ_T :: T  = 0.5   # Exponent for effect of unstable buoyancy forcing on wind mixing of T
    Cmb_U :: T  = 1/3   # Exponent for the effect of wind on convective mixing of U
    Cmb_T :: T  = 1/3   # Exponent for effect of wind on convective mixing of T

    Cd_U  :: T  = 0.5   # Wind mixing regime threshold for momentum
    Cd_T  :: T  = 2.5   # Wind mixing regime threshold for tracers

    Cb_U  :: T  = 0.599 # Buoyancy flux parameter for convective turbulence
    Cb_T  :: T  = 1.36  # Buoyancy flux parameter for convective turbulence
    Cτb_U :: T  = (Cτ / Cb_U)^(1/Cmb_U) * (1 + Cunst*Cd_U)^(Cmτ_U/Cmb_U) - Cd_U  # Wind stress parameter for convective turbulence
    Cτb_T :: T  = (Cτ / Cb_T)^(1/Cmb_T) * (1 + Cunst*Cd_T)^(Cmτ_T/Cmb_T) - Cd_T  # Wind stress parameter for convective turbulence

    CRi   :: T  = 0.3   # Critical bulk Richardson number
    CKE   :: T  = 4.32  # Unresolved turbulence parameter
    CKE₀  :: T  = 1e-11 # Minimum unresolved turbulence kinetic energy

    KU₀   :: T  = 1e-6  # Interior viscosity for velocity
    KT₀   :: T  = 1e-7  # Interior diffusivity for temperature
    KS₀   :: T  = 1e-9  # Interior diffusivity for salinity
end

# Shape functions.
# 'd' is a non-dimensional depth coordinate.
default_NL_shape(d) = ifelse(0<d<1, d*(1-d)^2, -zero(d))
const default_K_shape = default_NL_shape

mutable struct State{T} <: FieldVector{6, T}
    Fu :: T
    Fv :: T
    Fθ :: T
    Fs :: T
    Fb :: T
    h  :: T
end

State(T=Float64) = State{T}(0, 0, 0, 0, 0, 0)


mutable struct Model{S, G, T, U, B} <: AbstractModel{S, G, T}
    clock       :: Clock{T}
    grid        :: G
    timestepper :: S
    solution    :: U
    bcs         :: B
    parameters  :: Parameters{T}
    constants   :: Constants{T}
    state       :: State{T}
end


function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
      parameters = Parameters(),
         stepper = :ForwardEuler
              )

     K = (U=KU, V=KV, T=KT, S=KS)
     R = (U=RU, V=RV, T=RT, S=RS)
    eq = Equation(R=R, K=K, update=update_state!)

    bcs = (
        U = DefaultBoundaryConditions(eltype(grid)),
        V = DefaultBoundaryConditions(eltype(grid)),
        T = DefaultBoundaryConditions(eltype(grid)),
        S = DefaultBoundaryConditions(eltype(grid))
    )

    solution = Solution(
        CellField(grid),
        CellField(grid),
        CellField(grid),
        CellField(grid)
    )

    lhs = OceanTurb.build_lhs(solution)
    timestepper = Timestepper(stepper, eq, solution, lhs)
    clock = Clock()
    state = State()

    return Model(clock, grid, timestepper, solution, bcs, parameters, constants, state)
end

# Note: we use 'm' to refer to 'model' in function definitions below.

@inline Fb(g, α, β, Fθ, Fs) = g * (α*Fθ - β*Fs)

@propagate_inbounds d(m, i) = ifelse(m.state.h>0, -m.grid.zf[i]/m.state.h, -zero(m.state.h))

"Return the buoyancy gradient at face point i."
@propagate_inbounds ∂B∂z(T, S, g, α, β, i) = g * (α*∂z(T, i) - β*∂z(S, i))

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
    m.state.Fb = Fb(m.constants.g, m.constants.α, m.constants.β, m.state.Fθ, m.state.Fs)
    m.state.h  = mixing_depth(m)
    return nothing
end

#
# Diagnosis of mixing depth "h"
#

"Returns the surface_layer_average for mixing depth h = -zf[i]."
@propagate_inbounds function surface_layer_average(c, CSL::T, i) where T
    if i > c.grid.N # Return surface value
        return onface(c, c.grid.N+1)
    else
        iε = length(c)+1 - CSL*(length(c)+1 - i) # (fractional) face "index" of the surface layer
        face = ceil(Int, iε)  # next cell face above the fractional depth
        frac = face - iε # fraction of lowermost cell in the surface layer.
        surface_layer_integral = zero(T)

        # Contribution of fractional cell to total integral
        if face > 1
            surface_layer_integral += frac * Δf(c, face-1) * c[face-1]
        else
            face = 1
        end

        # Add cells above face, if there are any.
        for j = face:length(c)
            surface_layer_integral += Δf(c, j) * c[j]
        end

        h = -c.grid.zf[i] # depth

        return surface_layer_integral / (CSL*h)
    end
end

"""
Return Δc(hᵢ), the difference between the surface-layer average of c and its value at depth hᵢ, where
i is a face index.
"""
@propagate_inbounds Δ(c, CSL, i) = surface_layer_average(c, CSL, i) - onface(c, i)

"Returns the parameterization for unresolved KE at face point i."
@inline function unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, β)
    return CKE * h^(4/3) * sqrt(max(0, Bz)) * max(0, Fb)^(1/3) + CKE₀
end

"""
    bulk_richardson_number(model, i)

Returns the bulk Richardson number of `model` at face `i`.
"""
@propagate_inbounds function bulk_richardson_number(
            U, V, T, S, Fb::TT, CKE::TT, CKE₀::TT, CSL::TT,
            g::TT, α::TT, β::TT, i) where TT

    h = -U.grid.zf[i]
    # (h - hε) * ΔB
    h⁺ΔB = h * (one(TT) - CSL/2) * g * (α*Δ(T, CSL, i) - β*Δ(S, CSL, i))

    KE = (Δ(U, CSL, i)^2 + Δ(V, CSL, i)^2
              + unresolved_kinetic_energy(h, ∂B∂z(T, S, g, α, β, i), Fb, CKE, CKE₀, g, α, β))

    if KE == 0 && h⁺ΔB == 0 # Alistar Adcroft's theorem
        return -zero(TT)
    else
        return h⁺ΔB / KE
    end
end

@propagate_inbounds bulk_richardson_number(m, i) = bulk_richardson_number(
    m.solution.U, m.solution.V, m.solution.T, m.solution.S,
    m.state.Fb, m.parameters.CKE, m.parameters.CKE₀, m.parameters.CSL, m.constants.g,
    m.constants.α, m.constants.β, i)

"""
    mixing_depth(model)

Calculate the mixing depth 'h' for `model`.
"""
function mixing_depth(m)
    ih₁ = m.grid.N + 1 # start at top.
    @inbounds Ri₁ = bulk_richardson_number(m, ih₁) # should be 0.

    # Descend through grid until Ri rises above critical value
    while ih₁ > 1 && Ri₁ < m.parameters.CRi
        ih₁ -= 1 # descend
        @inbounds Ri₁ = bulk_richardson_number(m, ih₁)
    end

    # Edge cases:
    # 1. Mixing depth is at the top of the domain (z=0):
    if ih₁ == m.grid.N + 1
        @inbounds z★ = m.grid.zf[ih₁]

    # 2. Mixing depth is whole domain because Ri is always less than CRi:
    elseif ih₁ == 1 && Ri₁ < m.parameters.CRi
        @inbounds z★ = m.grid.zf[ih₁]

    # 3. Ri is infinite somewhere inside the domain.
    elseif !isfinite(Ri₁)
        @inbounds z★ = m.grid.zc[ih₁]

    # Main case: mixing depth is in the interior.
    else # Ri₁ > CRi
        ΔRi = bulk_richardson_number(m, ih₁+1) - Ri₁ # <0 linearly interpolate to find h.
        # x = x₀ + Δx * (y-y₀) / Δy
        @inbounds z★ = m.grid.zf[ih₁] + Δf(m.grid, ih₁) * (m.parameters.CRi - Ri₁) / ΔRi
    end

    -z★ < 0 && @warn "mixing depth $(-z★) is negative"

    return -z★ # "depth" is negative height.
end

#
# Vertical velocity scale
#

"Return true if the boundary layer is unstable and convecting."
@inline isunstable(model) = model.state.Fb > 0

"Return true if the boundary layer is forced."
@inline isforced(model) = model.state.Fu != 0 || model.state.Fv != 0 || model.state.Fb != 0

"Return the turbuent velocity scale associated with wind stress."
@inline ωτ(Fu, Fv) = (Fu^2 + Fv^2)^(1/4)
@inline ωτ(m::AbstractModel) = ωτ(m.state.Fu, m.state.Fv)

"Return the turbuent velocity scale associated with convection."
@inline ωb(Fb, h) = abs(h * Fb)^(1/3)
@inline ωb(m::AbstractModel) = ωb(m.state.Fb, m.state.h)

"Return the vertical velocity scale at depth d for a stable boundary layer."
@inline 𝒲_stable(Cτ, Cstab, Cn, ωτ, ωb, d) = Cτ * ωτ / (1 + Cstab * d * (ωb/ωτ)^3)^Cn

"Return the vertical velocity scale at scaled depth dϵ for an unstable boundary layer."
@inline function 𝒲_unstable(CSL, Cd, Cτ, Cunst, Cb, Cτb, Cmτ, Cmb, ωτ, ωb, d)
    dϵ = min(CSL, d)
    if dϵ * ωb^3 < Cd * ωτ^3
        return Cτ * ωτ * (1 + Cunst * dϵ * (ωb/ωτ)^3)^Cmτ
    else
        return Cb * ωb * (dϵ + Cτb * (ωτ/ωb)^3)^Cmb
    end
end

@propagate_inbounds function 𝒲_unstable_U(m, i)
    return 𝒲_unstable(m.parameters.CSL, m.parameters.Cd_U,
                            m.parameters.Cτ, m.parameters.Cunst,
                            m.parameters.Cb_U, m.parameters.Cτb_U,
                            m.parameters.Cmτ_U, m.parameters.Cmb_U,
                            ωτ(m), ωb(m), d(m, i)
                            )
end

@propagate_inbounds function 𝒲_unstable_T(m, i)
    return 𝒲_unstable(m.parameters.CSL, m.parameters.Cd_T,
                            m.parameters.Cτ, m.parameters.Cunst,
                            m.parameters.Cb_T, m.parameters.Cτb_T,
                            m.parameters.Cmτ_T, m.parameters.Cmb_T,
                            ωτ(m), ωb(m), d(m, i)
                            )
end

@propagate_inbounds function 𝒲_stable(m, i)
    return 𝒲_stable(m.parameters.Cτ, m.parameters.Cstab, m.parameters.Cn,
                          ωτ(m), ωb(m), d(m, i)
                          )
end

"Return the turbulent velocity scale for momentum at face point i."
@propagate_inbounds function 𝒲_U(m::AbstractModel{TS, G, T}, i) where {TS, G, T}
    if !isforced(m)
        return -zero(T)
    elseif isunstable(m)
        return 𝒲_unstable_U(m, i)
    else
        return 𝒲_stable(m, i)
    end
end

"Return the turbulent velocity scale for tracers at face point i."
@propagate_inbounds function 𝒲_T(m::AbstractModel{TS, G, T}, i) where {TS, G, T}
    if !isforced(m)
        return -zero(T)
    elseif isunstable(m)
        return 𝒲_unstable_T(m, i)
    else
        return 𝒲_stable(m, i)
    end
end

const 𝒲_V = 𝒲_U
const 𝒲_S = 𝒲_T

## ** The K-Profile-Parameterization **
K_KPP(h, 𝒲, d, shape=default_K_shape) = ifelse(0<d<1, max(zero(h), h*𝒲*shape(d)), -zero(h))

#
# Non-local flux
#

"""
    NL(CNL, flux, d, shape=default_shape)

Returns the nonlocal flux, N = CNL*flux*shape(d),
where `flux` is the flux of some quantity out of the surface,
`shape` is a shape function, and `d` is a non-dimensional depth coordinate
that increases from 0 at the surface to 1 at the bottom of the mixing layer.

Because flux is defined as pointing in the positive direction,
a positive surface flux implies negative surface flux divergence,
which implies a reduction to the quantity in question.
For example, positive heat flux out of the surface implies cooling.
"""
@inline NL(CNL, flux, d, shape=default_NL_shape) = CNL * flux * shape(d)

@inline function ∂NL∂z(CNL::T, Fϕ, dᵢ₊₁, dᵢ, Δf, m) where T
    if isunstable(m)
        return (NL(CNL, Fϕ, dᵢ₊₁) - NL(CNL, Fϕ, dᵢ)) / Δf
    else
        return -zero(T)
    end
end

@propagate_inbounds ∂NLT∂z(m, i) =
    ∂NL∂z(m.parameters.CNL, m.state.Fθ, d(m, i+1), d(m, i), Δf(m.grid, i), m)

@propagate_inbounds ∂NLS∂z(m, i) =
    ∂NL∂z(m.parameters.CNL, m.state.Fs, d(m, i+1), d(m, i), Δf(m.grid, i), m)

#
# Equation specification
#

# K_{U,V,T,S} is calculated at face points
@propagate_inbounds KU(m, i) = K_KPP(m.state.h, 𝒲_U(m, i), d(m, i)) + m.parameters.KU₀
@propagate_inbounds KT(m, i) = K_KPP(m.state.h, 𝒲_T(m, i), d(m, i)) + m.parameters.KT₀
@propagate_inbounds KS(m, i) = K_KPP(m.state.h, 𝒲_S(m, i), d(m, i)) + m.parameters.KS₀
const KV = KU

@propagate_inbounds RU(m, i) =   m.constants.f * m.solution.V[i]
@propagate_inbounds RV(m, i) = - m.constants.f * m.solution.U[i]
@propagate_inbounds RT(m, i) = - ∂NLT∂z(m, i)
@propagate_inbounds RS(m, i) = - ∂NLS∂z(m, i)

#####
##### Some utilities
#####

function nonlocal_salinity_flux!(flux, m)
    for i in interiorindices(flux)
        @inbounds flux[i] = NL(m.parameters.CNL, m.state.Fs, d(m, i))
    end
    return nothing
end

function nonlocal_temperature_flux!(flux, m)
    for i in interiorindices(flux)
        @inbounds flux[i] = NL(m.parameters.CNL, m.state.Fθ, d(m, i))
    end
    return nothing
end

function nonlocal_salinity_flux(model)
    flux = FaceField(model.grid)
    nonlocal_salinity_flux!(flux, model)
    return flux
end

function nonlocal_temperature_flux(model)
    flux = FaceField(model.grid)
    nonlocal_temperature_flux!(flux, model)
    return flux
end

end # module
