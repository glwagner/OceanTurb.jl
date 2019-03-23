module KPP

using
    OceanTurb,
    StaticArrays,
    LinearAlgebra

const nsol = 4
@specify_solution CellField U V T S

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
    CNL   :: T  # Non-local flux proportionality constant

    Cstab :: T  # Reduction of wind-driven diffusivity due to stable buoyancy flux
    Cunst :: T  # Reduction of wind-driven diffusivity due to stable buoyancy flux

    Cb_U  :: T  # Buoyancy flux viscosity proportionality for convective turbulence
    Cτ_U  :: T  # Wind stress viscosity proportionality for convective turbulence
    Cb_T  :: T  # Buoyancy flux diffusivity proportionality for convective turbulence
    Cτ_T  :: T  # Wind stress diffusivity proportionality for convective turbulence

    Cd_U  :: T  # Buoyancy flux diffusivity proportionality for convective turbulence
    Cd_T  :: T  # Wind stress diffusivity proportionality for convective turbulence

    CRi   :: T  # Critical bulk_richardson_number number
    CKE   :: T  # Unresolved turbulence parameter

    KU₀   :: T  # Interior diffusivity
    KT₀   :: T  # Interior diffusivity
    KS₀   :: T  # Interior diffusivity
end

function Parameters(T=Float64;
       Cε = 0.1,
       Cκ = 0.4,
      CNL = 6.33,
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
       K₀ = 1e-5, KU₀=K₀, KT₀=K₀, KS₀=K₀
     )

     Parameters{T}(Cε, Cκ, CNL, Cstab, Cunst,
                   Cb_U, Cτ_U, Cb_T, Cτ_T, Cd_U, Cd_T,
                   CRi, CKE, KU₀, KT₀, KS₀)
end

# Shape functions (these shoul become parameters eventually).
# 'd' is a non-dimensional depth coordinate.
default_shape_N(d) = d*(1-d)^2
default_shape_K(d) = d*(1-d)^2

struct Constants{T}
    g  :: T # Gravitiational acceleration
    cP :: T # Heat capacity of water
    ρ₀ :: T # Reference density
    α  :: T # Thermal expansion coefficient
    β  :: T # Haline expansion coefficient
    f  :: T # Coriolis parameter
end

function Constants(T=Float64; α=2.5e-4, β=8e-5, ρ₀=1035, cP=3992, f=0, g=9.81)
    Constants{T}(g, cP, ρ₀, α, β, f)
end

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

struct Model{TS, G, E, T} <: AbstractModel{TS, G, E, T}
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
    equation = Equation(calc_rhs_explicit!)
    timestepper = Timestepper(:ForwardEuler, solution)

    return Model(timestepper, grid, equation, solution, bcs, Clock(),
                    parameters, constants, State())
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

# K_{U,V,T,S} is calculated at face points
K_U(m, i) = K_KPP(m.state.h, w_scale_U(m, i), d(m, i)) + m.parameters.KU₀
K_T(m, i) = K_KPP(m.state.h, w_scale_T(m, i), d(m, i)) + m.parameters.KT₀
K_S(m, i) = K_KPP(m.state.h, w_scale_S(m, i), d(m, i)) + m.parameters.KS₀
const K_V = K_U

"Return the buoyancy gradient at face point i."
∂B∂z(T, S, g, α, β, i) = g * (α*∂z(T, i) - β*∂z(S, i))
∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α, m.constants.β, i)

#
# Diagnosis of mixing depth "h"
#

"Returns the surface_layer_average for mixing depth h = -zf[i]."
function surface_layer_average(c, Cε, i)
    iε = length(c)+1 - Cε*(length(c)+1 - i) # (fractional) face "index" of the surface layer
    face = ceil(Int, iε)  # the next cell face above the fractional depth
    frac = face - iε # the fraction of the lowest cell in the surface layer.

    # Example 1:

    #   length(c) = 9 (face_length = 10)
    #          Cε = 0.1
    #           i = 9
    #   => iε = 10 - 0.1*(1) = 9.9, face = 10, frac = 0.1.

    # Example 2:

    # length(c) = 99 (face_length = 100)
    #        Cε = 0.1
    #         i = 18
    #       => iε = 100 - 0.1*82 = 91.8, face = 92, frac = 0.2.

    # Contribution of fractional cell to total integral
    surface_layer_integral = frac > 0 ? frac * Δf(c, face-1) * c.data[face-1] : 0

    # Add cells above face, if there are any.
    for j = face:length(c)
      @inbounds surface_layer_integral += Δf(c, j) * c.data[j]
    end

    h = -c.grid.zf[i] # depth
    return surface_layer_integral / (Cε*h)
end

"""
Return Δc(hᵢ), the difference between the surface-layer average of c and its value at depth hᵢ, where
i is a face index.
"""
Δ(c::CellField, Cε, i) = surface_layer_average(c, Cε, i) - onface(c, i)

"Returns the parameterization for unresolved KE at face point i."
function unresolved_kinetic_energy(T, S, Bz, Fb, CKE, g, α, β, i)
    h = -T.grid.zf[i]
    return CKE * h^(4/3) * sqrt(max(0, Bz)) * max(0, Fb)^(1/3)
end

"""
    bulk_richardson_number(model, i)

Returns the bulk Richardson number of `model` at face `i`.
"""
function bulk_richardson_number(U, V, T, S, Fb, CKE, Cε, g, α, β, i)
    hΔB = -U.grid.zf[i] * ( g * (α*Δ(T, Cε, i) - β*Δ(S, Cε, i)) )
    Bz = ∂B∂z(T, S, g, α, β, i)
    uKE = unresolved_kinetic_energy(T, S, Bz, Fb, CKE, g, α, β, i)
    KE = Δ(U, Cε, i)^2 + Δ(V, Cε, i)^2 + uKE

    if KE == 0 && hΔB == 0 # Alistar Adcroft's theorem
        return 0
    else
        return hΔB / KE
    end
end

bulk_richardson_number(m, i) = bulk_richardson_number(
    m.solution.U, m.solution.V, m.solution.T, m.solution.S,
    m.state.Fb, m.parameters.CKE, m.parameters.Cε, m.constants.g,
    m.constants.α, m.constants.β, i)

"""
    mixing_depth(model)

Calculate the mixing depth 'h' for `model`.
"""
function mixing_depth(m)
    # Descend through grid until Ri rises above critical value
    Ri₁ = 0
    ih₁ = m.grid.N + 1 # start at top
    while ih₁ > 2 && Ri₁ < m.parameters.CRi
        ih₁ -= 1 # descend
        Ri₁ = bulk_richardson_number(m, ih₁)
    end

    # Here, ih₁ >= 2.

    if !isfinite(Ri₁)         # Ri is infinite:
        z★ = m.grid.zf[ih₁+1] # "mixing depth" is just above where Ri = inf.

    elseif Ri₁ < m.parameters.CRi # We descended to ih₁=2 and Ri is still too low:
        z★ = m.grid.zf[1]         # mixing depth extends to bottom of grid.

    else                                             # We have descended below critical Ri:
        if ih₁ == m.grid.N
            z★ = m.grid.zf[m.grid.N]
        else
            ΔRi = bulk_richardson_number(m, ih₁+1) - Ri₁ # linearly interpolate to find h.
            # x = x₀ + Δx * (y-y₀) / Δy
            z★ = m.grid.zf[ih₁] + Δf(m.grid, ih₁) * (m.parameters.CRi - Ri₁) / ΔRi
        end
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

const nU = 1/4
const nT = 1/2

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
    nonlocal_flux(flux, d, shape=default_shape)

Returns the nonlocal flux, N = flux*shape(d),
where `flux` is the flux of some quantity out of the surface,
`shape` is a shape function, and `d` is a non-dimensional depth coordinate
that increases from 0 at the surface to 1 at the bottom of the mixing layer.

Because flux is defined as pointing in the positive direction,
a positive surface flux implies negative surface flux divergence,
which implies a reduction to the quantity in question.
For example, positive heat flux out of the surface implies cooling.
"""
function nonlocal_flux(flux, d, shape=default_shape_N)
    if 0 < d < 1
        return flux*shape(d) # not minus sign due to flux convention
    else
        return 0
    end
end

const N = nonlocal_flux

∂NT∂z(m, i) = ( N(m.state.Fθ, d(m, i+1)) - N(m.state.Fθ, d(m, i)) ) / Δf(m.grid, i)
∂NS∂z(m, i) = ( N(m.state.Fs, d(m, i+1)) - N(m.state.Fs, d(m, i)) ) / Δf(m.grid, i)

#
# Local diffusive flux
#

const BC = BoundaryCondition

# ∇K∇c for c::CellField
K∂z(K, c, i) = K*∂z(c, i)
∇K∇c(Kᵢ₊₁, Kᵢ, c, i)              = ( K∂z(Kᵢ₊₁, c, i+1) -    K∂z(Kᵢ, c, i)      ) /    Δf(c, i)
∇K∇c_top(Kᵢ, c, top_flux)         = (     -top_flux     - K∂z(Kᵢ, c, length(c)) ) / Δf(c, length(c))
∇K∇c_bottom(Kᵢ₊₁, c, bottom_flux) = (  K∂z(Kᵢ₊₁, c, 2)  +     bottom_flux       ) /    Δf(c, 1)

## Top and bottom flux estimates for constant (Dirichlet) boundary conditions
bottom_flux(K, c, c_bndry, Δf) = -2K*( bottom(c) - c_bndry ) / bottom(Δf) # -K*∂c/∂z at the bottom
top_flux(K, c, c_bndry, Δf)    = -2K*(  c_bndry  -  top(c) ) /   top(Δf)  # -K*∂c/∂z at the top

∇K∇c_top(Kᵢ, c, bc::BC{<:Flux}, model) = ∇K∇c_top(Kᵢ, c, get_bc(bc, model))
∇K∇c_bottom(Kᵢ₊₁, Kᵢ, c, bc::BC{<:Flux}, model) = ∇K∇c_bottom(Kᵢ₊₁, c, getbc(model, bc))
∇K∇c_bottom(Kᵢ₊₁, Kᵢ, c, bc::BC{<:Gradient}, model) = ∇K∇c_bottom(Kᵢ₊₁, c, -Kᵢ*getbc(model, bc))

function ∇K∇c_bottom(Kᵢ₊₁, Kᵢ, c, bc::BC{<:Value}, model)
    flux = bottom_flux(Kᵢ, c, getbc(model, bc), Δf(model.grid, 1))
    return ∇K∇c_bottom(Kᵢ₊₁, c, flux)
end

#
# Equation entry
#

function calc_rhs_explicit!(∂t, m)

    # Preliminaries
    U, V, T, S = m.solution
    update_state!(m)

    for i in interior(U)
        @inbounds begin
            ∂t.U[i] = ∇K∇c(K_U(m, i+1), K_U(m, i), U, i) + m.constants.f*V[i]
            ∂t.V[i] = ∇K∇c(K_V(m, i+1), K_V(m, i), V, i) - m.constants.f*U[i]
            ∂t.T[i] = ∇K∇c(K_T(m, i+1), K_T(m, i), T, i) - ∂NT∂z(m, i)
            ∂t.S[i] = ∇K∇c(K_S(m, i+1), K_S(m, i), S, i) - ∂NS∂z(m, i)
        end
    end

    # Flux into the top (the only boundary condition allowed)
    i = m.grid.N
    ∂t.U[i] = ∇K∇c_top(K_U(m, i), U, m.state.Fu) + m.constants.f*V[i]
    ∂t.V[i] = ∇K∇c_top(K_V(m, i), V, m.state.Fv) - m.constants.f*U[i]
    ∂t.T[i] = ∇K∇c_top(K_T(m, i), T, m.state.Fθ) - ∂NT∂z(m, i)
    ∂t.S[i] = ∇K∇c_top(K_S(m, i), S, m.state.Fs) - ∂NS∂z(m, i)

    # Bottom
    i = 1
    ∂t.U[i] = ∇K∇c_bottom(K_U(m, i+1), K_U(m, i), U, m.bcs.U.bottom, m) + m.constants.f*V[i]
    ∂t.V[i] = ∇K∇c_bottom(K_V(m, i+1), K_V(m, i), V, m.bcs.V.bottom, m) - m.constants.f*U[i]
    ∂t.T[i] = ∇K∇c_bottom(K_T(m, i+1), K_T(m, i), T, m.bcs.T.bottom, m) - ∂NT∂z(m, i)
    ∂t.S[i] = ∇K∇c_bottom(K_S(m, i+1), K_S(m, i), S, m.bcs.S.bottom, m) - ∂NS∂z(m, i)

    return nothing
end

end # module
