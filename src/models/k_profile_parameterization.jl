module KPP

using
    OceanTurb,
    LinearAlgebra

import OceanTurb: Constants

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

Fb(g, α, β, Fθ, Fs) = g * (α*Fθ - β*Fs)

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

mutable struct Model{TS, G, T} <: AbstractModel{TS, G, T}
    @add_standard_model_fields
    parameters :: Parameters{T}
    constants  :: Constants{T}
    state      :: State{T}
end

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
      parameters = Parameters(),
         stepper = :ForwardEuler,
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
K_KPP(h, w_scale, d, shape=default_shape_K) = 0 < d < 1 ? max(0, h * w_scale * shape(d)) : 0
d(m, i) = -m.grid.zf[i] / m.state.h

"Return the buoyancy gradient at face point i."
∂B∂z(T, S, g, α, β, i) = g * (α*∂z(T, i) - β*∂z(S, i))
∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α, m.constants.β, i)

#
# Diagnosis of mixing depth "h"
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
    ih₁ = m.grid.N + 1 # start at top.
    Ri₁ = bulk_richardson_number(m, ih₁) # should be 0.

    # Descend through grid until Ri rises above critical value
    while ih₁ > 1 && Ri₁ < m.parameters.CRi
        ih₁ -= 1 # descend
        Ri₁ = bulk_richardson_number(m, ih₁)
    end

    # Edge cases:
    # 1. Mixing depth is 0:
    if ih₁ == m.grid.N + 1
        z★ = m.grid.zf[ih₁]

    # 2. Mixing depth is whole domain because Ri is always less than CRi:
    elseif ih₁ == 1 && Ri₁ < m.parameters.CRi
        z★ = m.grid.zf[ih₁]

    # 3. Ri is infinite somewhere inside the domain.
    elseif !isfinite(Ri₁)
        z★ = m.grid.zc[ih₁]

    # Main case: mixing depth is in the interior.
    else # Ri₁ > CRi
        ΔRi = bulk_richardson_number(m, ih₁+1) - Ri₁ # <0 linearly interpolate to find h.
        # x = x₀ + Δx * (y-y₀) / Δy
        z★ = m.grid.zf[ih₁] + Δf(m.grid, ih₁) * (m.parameters.CRi - Ri₁) / ΔRi
    end

    -z★ < 0 && @warn "mixing depth $(-z★) is negative"

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
w_scale_stable(Cτ, Cstab, Cn, ωτ, ωb, d) = Cτ * ωτ / (1 + Cstab * d * (ωb/ωτ)^3)^Cn

"Return the vertical velocity scale at scaled depth dϵ for an unstable boundary layer."
function w_scale_unstable(CSL, Cd, Cτ, Cunst, Cb, Cτb, Cmτ, Cmb, ωτ, ωb, d)
    dϵ = min(CSL, d)
    if dϵ < Cd * (ωτ/ωb)^3
        return Cτ * ωτ * (1 + Cunst * dϵ * (ωb/ωτ)^3)^Cmτ
    else
        return Cb * ωb * (dϵ + Cτb * (ωτ/ωb)^3)^Cmb
    end
end

function w_scale_unstable_U(m, i)
    return w_scale_unstable(m.parameters.CSL, m.parameters.Cd_U,
                            m.parameters.Cτ, m.parameters.Cunst,
                            m.parameters.Cb_U, m.parameters.Cτb_U,
                            m.parameters.Cmτ_U, m.parameters.Cmb_U,
                            ωτ(m), ωb(m), d(m, i)
                            )
end

function w_scale_unstable_T(m, i)
    return w_scale_unstable(m.parameters.CSL, m.parameters.Cd_T,
                            m.parameters.Cτ, m.parameters.Cunst,
                            m.parameters.Cb_T, m.parameters.Cτb_T,
                            m.parameters.Cmτ_T, m.parameters.Cmb_T,
                            ωτ(m), ωb(m), d(m, i)
                            )
end

function w_scale_stable(m, i)
    return w_scale_stable(m.parameters.Cτ, m.parameters.Cstab, m.parameters.Cn,
                          ωτ(m), ωb(m), d(m, i)
                          )
end

"Return the vertical velocity scale for momentum at face point i."
function w_scale_U(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return w_scale_unstable_U(m, i)
    else
        return w_scale_stable(m, i)
    end
end


"Return the vertical velocity scale for tracers at face point i."
function w_scale_T(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return w_scale_unstable_T(m, i)
    else
        return w_scale_stable(m, i)
    end
end

const w_scale_V = w_scale_U
const w_scale_S = w_scale_T

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
N(CNL, flux, d, shape=default_shape_N) = CNL * flux * shape(d)

function ∂N∂z(CNL, Fϕ, d, Δf, m)
    if isunstable(m)
        return (N(CNL, Fϕ, d) - N(CNL, Fϕ, d)) / Δf
    else
        return 0
    end
end

∂NT∂z(m, i) = @inbounds ∂N∂z(m.parameters.CNL, m.state.Fθ, d(m, i), Δf(m.grid, i), m)
∂NS∂z(m, i) = @inbounds ∂N∂z(m.parameters.CNL, m.state.Fs, d(m, i), Δf(m.grid, i), m)

#
# Equation specification
#

# K_{U,V,T,S} is calculated at face points
KU(m, i) = K_KPP(m.state.h, w_scale_U(m, i), d(m, i)) + m.parameters.KU₀
KT(m, i) = K_KPP(m.state.h, w_scale_T(m, i), d(m, i)) + m.parameters.KT₀
KS(m, i) = K_KPP(m.state.h, w_scale_S(m, i), d(m, i)) + m.parameters.KS₀
const KV = KU

@inline RU(f, V, i) = @inbounds  f*V[i]
@inline RV(f, U, i) = @inbounds -f*U[i]

@inline RU(m, i) = RU(m.constants.f, m.solution.V, i)
@inline RV(m, i) = RV(m.constants.f, m.solution.U, i)
@inline RT(m, i) = -∂NT∂z(m, i)
@inline RS(m, i) = -∂NS∂z(m, i)

end # module
