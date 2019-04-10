module KPP

using
    OceanTurb,
    LinearAlgebra

import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top, Constants

const nU = 1/4 # exponent for momentum turbulent velocity scale
const nT = 1/2 # exponent for tracer turbulent velocity scale

const nsol = 4
@solution U V T S

"""
    Parameters(; kwargs...)

Construct KPP parameters.

    Args
    ====
    Cε : Surface layer fraction
    etc.
"""
struct Parameters{T} <: AbstractParameters
    Cε    :: T  # Surface layer fraction
    Cκ    :: T  # Von Karman constant
    CN    :: T  # Non-local flux proportionality constant

    Cstab :: T  # Stable buoyancy flux parameter for wind-driven turbulence
    Cunst :: T  # Unstable buoyancy flux parameter for wind-driven turbulence

    Cb_U  :: T  # Buoyancy flux parameter for convective turbulence
    Cτ_U  :: T  # Wind stress parameter for convective turbulence
    Cb_T  :: T  # Buoyancy flux parameter for convective turbulence
    Cτ_T  :: T  # Wind stress parameter for convective turbulence

    Cd_U  :: T  # Wind mixing regime threshold for momentum
    Cd_T  :: T  # Wind mixing regime threshold for tracers

    CRi   :: T  # Critical bulk Richardson number
    CKE   :: T  # Unresolved turbulence parameter
    CKE₀  :: T  # Minimum unresolved turbulence kinetic energy

    KU₀   :: T  # Interior viscosity for velocity
    KT₀   :: T  # Interior diffusivity for temperature
    KS₀   :: T  # Interior diffusivity for salinity
end

function Parameters(T=Float64;
       Cε = 0.1,
       Cκ = 0.4,
       CN = 6.33,
    Cstab = 2.0,
    Cunst = 6.4,
     Cb_U = 0.599,
     Cτ_U = 0.135,
     Cb_T = 1.36,
     Cτ_T = -1.85,
     Cd_U = 0.5,
     Cd_T = 2.5,
      CRi = 4.32,
      CKE = 0.3,
     CKE₀ = 1e-11,
       K₀ = 1e-5, KU₀=K₀, KT₀=K₀, KS₀=K₀
     )

     Parameters{T}(Cε, Cκ, CN, Cstab, Cunst,
                   Cb_U, Cτ_U, Cb_T, Cτ_T, Cd_U, Cd_T,
                   CRi, CKE, CKE₀, KU₀, KT₀, KS₀)
end

# Shape functions (these shoul become parameters eventually).
# 'd' is a non-dimensional depth coordinate.
default_shape_N(d) = d*(1-d)^2
default_shape_K(d) = d*(1-d)^2

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

struct Model{TS, G, T} <: AbstractModel{TS, G, T}
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
function K_KPP(h, w_scale, d, shape=default_shape_K)
    if 0 < d < 1
        return max(0, h * w_scale * shape(d))
    else
        return 0.0
    end
end

d(m, i) = -m.grid.zf[i] / m.state.h

"Return the buoyancy gradient at face point i."
∂B∂z(T, S, g, α, β, i) = g * (α*∂z(T, i) - β*∂z(S, i))
∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α, m.constants.β, i)

#
# Diagnosis of mixing depth "h"
#

"Returns the surface_layer_average for mixing depth h = -zf[i]."
function surface_layer_average(c, Cε, i)
    if i > c.grid.N # Return surface value
        return onface(c, c.grid.N+1)
    else
        iε = length(c)+1 - Cε*(length(c)+1 - i) # (fractional) face "index" of the surface layer
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

        return surface_layer_integral / (Cε*h)
    end
end

"""
Return Δc(hᵢ), the difference between the surface-layer average of c and its value at depth hᵢ, where
i is a face index.
"""
Δ(c, Cε, i) = surface_layer_average(c, Cε, i) - onface(c, i)

"Returns the parameterization for unresolved KE at face point i."
function unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, β, i)
    return CKE * h^(4/3) * sqrt(max(0, Bz)) * max(0, Fb)^(1/3) + CKE₀
end

"""
    bulk_richardson_number(model, i)

Returns the bulk Richardson number of `model` at face `i`.
"""
function bulk_richardson_number(U, V, T, S, Fb, CKE, CKE₀, Cε, g, α, β, i)
    h = -U.grid.zf[i]
    # (h - hε) * ΔB
    h⁺ΔB = h * (1 - 0.5Cε) * g * (α*Δ(T, Cε, i) - β*Δ(S, Cε, i))

    Bz = ∂B∂z(T, S, g, α, β, i)
    unresolved_KE = unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, β, i)
    KE = Δ(U, Cε, i)^2 + Δ(V, Cε, i)^2 + unresolved_KE

    if KE == 0 && h⁺ΔB == 0 # Alistar Adcroft's theorem
        return 0
    else
        return h⁺ΔB / KE
    end
end

bulk_richardson_number(m, i) = bulk_richardson_number(
    m.solution.U, m.solution.V, m.solution.T, m.solution.S,
    m.state.Fb, m.parameters.CKE, m.parameters.CKE₀, m.parameters.Cε, m.constants.g,
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
isforced(model) = ωτ(model) > 0 || ωb(model) > 0

"Return the turbuent velocity scale associated with wind stress."
ωτ(Fu, Fv) = (Fu^2 + Fv^2)^(1/4)
ωτ(m::Model) = ωτ(m.state.Fu, m.state.Fv)

"Return the turbuent velocity scale associated with convection."
ωb(Fb, h) = abs(h * Fb)^(1/3)
ωb(m::Model) = ωb(m.state.Fb, m.state.h)

"Return truncated, non-dimensional depth coordinate."
dϵ(m::Model, d) = min(m.parameters.Cε, d)

"Return the vertical velocity scale at depth d for a stable boundary layer."
w_scale_stable(Cκ, Cstab, ωτ, ωb, d) = Cκ * ωτ / (1 + Cstab * d * (ωb/ωτ)^3)

"Return the vertical velocity scale at scaled depth dϵ for an unstable boundary layer."
function w_scale_unstable(Cd, Cκ, Cunst, Cb, Cτ, ωτ, ωb, dϵ, n)
    if dϵ < Cd * (ωτ/ωb)^3
        return Cκ * ωτ * ( 1 + Cunst * dϵ * (ωb/ωτ)^3 )^n
    else
        return Cb * ωb * ( dϵ + Cτ * (ωτ/ωb)^3 )^(1/3)
    end
end


"Return the vertical velocity scale for momentum at face point i."
function w_scale_U(m, i)
    if !isforced(m)
        return 0
    else
        if isunstable(m)
            return w_scale_unstable(m.parameters.Cd_U, m.parameters.Cκ, m.parameters.Cunst,
                                    m.parameters.Cb_U, m.parameters.Cτ_U,
                                    ωτ(m), ωb(m), min(m.parameters.Cε, d(m, i)), nU)
        else
            return w_scale_stable(m.parameters.Cκ, m.parameters.Cstab, ωτ(m), ωb(m), d(m, i))
        end
    end
end

"Return the vertical velocity scale for tracers at face point i."
function w_scale_T(m, i)
    if !isforced(m)
        return 0
    else
        if isunstable(m)
            return w_scale_unstable(m.parameters.Cd_T, m.parameters.Cκ, m.parameters.Cunst,
                                    m.parameters.Cb_T, m.parameters.Cτ_T,
                                    ωτ(m), ωb(m), min(m.parameters.Cε, d(m, i)), nT)
        else
            return w_scale_stable(m.parameters.Cκ, m.parameters.Cstab, ωτ(m), ωb(m), d(m, i))
        end
    end
end

const w_scale_V = w_scale_U
const w_scale_S = w_scale_T

#
# Non-local flux
#

"""
    N(CN, flux, d, shape=default_shape)

Returns the nonlocal flux, N = CN*flux*shape(d),
where `flux` is the flux of some quantity out of the surface,
`shape` is a shape function, and `d` is a non-dimensional depth coordinate
that increases from 0 at the surface to 1 at the bottom of the mixing layer.

Because flux is defined as pointing in the positive direction,
a positive surface flux implies negative surface flux divergence,
which implies a reduction to the quantity in question.
For example, positive heat flux out of the surface implies cooling.
"""
N(CN, flux, d, shape=default_shape_N) = 0 < d < 1 ? -CN*flux*shape(d) : 0

function ∂N∂z(CN, Fϕ, m, i)
    if isunstable(m)
        return (N(CN, Fϕ, d(m, i+1)) - N(CN, Fϕ, d(m, i))) / Δf(m.grid, i)
    else
        return 0
    end
end

∂NT∂z(m, i) = ∂N∂z(m.parameters.CN, m.state.Fθ, m, i)
∂NS∂z(m, i) = ∂N∂z(m.parameters.CN, m.state.Fs, m, i)

#
# Equation specification
#

# K_{U,V,T,S} is calculated at face points
KU(m, i) = K_KPP(m.state.h, w_scale_U(m, i), d(m, i)) + m.parameters.KU₀
KT(m, i) = K_KPP(m.state.h, w_scale_T(m, i), d(m, i)) + m.parameters.KT₀
KS(m, i) = K_KPP(m.state.h, w_scale_S(m, i), d(m, i)) + m.parameters.KS₀
const KV = KU

RU(m, i) =   m.constants.f * m.solution.V[i]
RV(m, i) = - m.constants.f * m.solution.U[i]
RT(m, i) = - ∂NT∂z(m, i)
RS(m, i) = - ∂NS∂z(m, i)

end # module
