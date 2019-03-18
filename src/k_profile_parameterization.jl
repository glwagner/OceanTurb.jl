module KPP

using
    Reexport,
    StaticArrays,
    LinearAlgebra

@reexport using OceanTurb

export
    Parameters,
    Constants,
    Model

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
    CNL   :: T  # Non-local flux proportionality constant
    Cκ    :: T  # Von Karman constant
    Cτb   :: T  # Reduction of wind-driven diffusivity due to stable buoyancy flux
    Cb_U  :: T  # Buoyancy flux viscosity proportionality for convective turbulence
    Cτ_U  :: T  # Wind stress viscosity proportionality for convective turbulence
    Cb_T  :: T  # Buoyancy flux diffusivity proportionality for convective turbulence
    Cτ_T  :: T  # Wind stress diffusivity proportionality for convective turbulence
    Cb⁺_U :: T  # Upper layer buoyancy flux viscosity proportionality for convective turbulence
    Cτ⁺_U :: T  # Upper layer wind stress viscosity proportionality for convective turbulence
    Cb⁺_T :: T  # Upper layer buoyancy flux diffusivity proportionality for convective turbulence
    Cτ⁺_T :: T  # Upper layer wind stress diffusivity proportionality for convective turbulence
    CRi   :: T  # Critical Richardson number
    CKE   :: T  # Unresolved turbulence parameter
    KU₀   :: T  # Interior diffusivity
    KT₀   :: T  # Interior diffusivity
    KS₀   :: T  # Interior diffusivity
end

function Parameters(T=Float64;
       Cε = 0.1,
      CNL = 6.33,
       Cκ = 0.4,
      Cτb = 2.0,
     Cb_U = 0.215,
     Cτ_U = 2.53,
     Cb_T = 0.0806,
     Cτ_T = -1.85,
    Cb⁺_U = 0.0256,
    Cτ⁺_U = 0.164,
    Cb⁺_T = 0.160,
    Cτ⁺_T = 1.02,
      CRi = 4.32,
      CKE = 0.3,
       K₀ = 1e-5, KU₀=K₀, KT₀=K₀, KS₀=K₀
     )

  Parameters{T}(Cε, CNL, Cκ, Cτb,
                  Cb_U, Cτ_U, Cb_T, Cτ_T,
                  Cb⁺_U, Cτ⁺_U, Cb⁺_T, Cτ⁺_T,
                  CRi, CKE, KU₀, KT₀, KS₀)
end


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

struct Model{TS, G, E, T, TP, TC, A} <: AbstractModel{TS, G, E, T}
    @add_standard_model_fields
    parameters::Parameters{TP}
    constants::Constants{TC}
    δRi::FaceField{A, G}
end

function Model(;
               N = 10,
               L = 1.0,
         stepper = :ForwardEuler,
             bcs = BoundaryConditions((ZeroFlux() for i=1:nsol)...),
       constants = Constants(),
      parameters = Parameters()
)

    grid = UniformGrid(N, L)
    solution = Solution((CellField(grid) for i=1:nsol)...)
    equation = Equation(calc_rhs_explicit!)
    timestepper = Timestepper(:ForwardEuler, solution)
    δRi = FaceField(grid)

    return Model(timestepper, grid, equation, solution, bcs, Clock(), parameters, constants, δRi)
end

#
# Equation specification
#
# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

 α(m) = m.constants.α
 β(m) = m.constants.β
T₀(m) = m.constants.T₀
S₀(m) = m.constants.S₀
ρ₀(m) = m.constants.ρ₀
 g(m) = m.constants.g

# Aliases for surface fluxes
Fu(m) = m.bcs.U.top.condition(m)
Fv(m) = m.bcs.V.top.condition(m)
Fθ(m) = m.bcs.T.top.condition(m)
Fs(m) = m.bcs.S.top.condition(m)
Fb(m) = g(m)/ρ₀(m) * (α(m)*Fθ(m) - β(m)*Fs(m))

const buoyancy_flux = Fb
∂B∂z(g, α, β, T, S, i) = g * (α*∂z(T, i) - β*∂z(S, i))

# Shape function. 'd' is a non-dimensional depth coordinate.
default_shape_N(d) = d*(1-d)^2
default_shape_K(d) = d*(1-d)^2

## Diffusivity

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
Δ(c, Cε, i) = surface_layer_average(c, Cε, i) - avz(c, i)

"Returns the parameterization for unresolved KE at face point i."
function unresolved_kinetic_energy(CKE, Fb, g, α, β, T, S, i)
  h = -T.grid.zf[i]
  return CKE * h^(4/3) * sqrt(max(0, ∂B∂z(g, α, β, T, S, i))) * max(0, Fb^(1/3))
end

"""
    Richardson(model, i)

Returns the Richardson number of `model`, including an estimate
for unresolved kinetic energy, at face point i.
"""
function Richardson(CKE, Cε, Fb, g, ρ₀, α, β, U, V, T, S, i)
  h = -T.grid.zf[i]
   ΔB = g/ρ₀ * (α*Δ(T, Cε, i) - β*Δ(S, Cε, i))
  ΔU² = Δ(U, Cε, i)^2 + Δ(V, Cε, i)^2
  return h*ΔB / (ΔU² + unresolved_kinetic_energy(CKE, Fb, g, α, β, T, S, i))
end

Richardson(m, i) = Richardson(m.parameters.CKE, m.parameters.Cε, buoyancy_flux(m),
                              m.constants.g, m.constants.ρ₀, m.constants.α, m.constants.β,
                              m.solution.U, m.solution.V, m.solution.T, m.solution.S, i)

"""
    mixing_depth(model)

Calculate 'h', the mixing depth, for the current state of `model`.
"""
function mixing_depth(m)
    # Descend through grid until Ri rises above critical value
    Ri₁ = 0
    #Ri₂ = 0
    ih₁ = m.grid.N
    while ih₁ > 1 && Ri₁ < m.parameters.CRi
        Ri₁ = Richardson(m, i)
        ih₁ -= 1
        #Ri₂ = Ri₁ < m.parameters.CRi ? Ri₁ : Ri₂
    end

    if Ri < m.parameters.CRi
        z★ = m.grid.zf[1] # mixing depth extends to bottom of grid
    else
        ih₂ = ih₁ + 1
        Δz = m.grid.zf[ih₂] - m.grid.zf[ih₁]
        Ri₂ = Richardson(m, ih₂)
        ΔRi = Ri₂ - Ri₁ # delta of Ri deviation

        # (x-x₀) = Δx * (y-y₀) / Δy
        z★ = m.grid.zf[ih₁] + Δz * (m.parameters.CRi - Ri₁) / ΔRi
    end

    return -z★
end

# Non-dimensional depth
face_depth(m, i) = -m.grid.zf[i] / mixing_depth(m)
cell_depth(m, i) = -m.grid.zc[i] / mixing_depth(m)

#
# Vertical velocity scale
#

ϖτ(Fu, Fv, ρ₀) = (Fu^2 + Fv^2)^(1/4)
ϖb(Fb, h) = abs(h * Fb)
dϵ(ϖτ, ϖb, d) = (ϖτ / ϖb)^3 * d

ϖτ(m::Model) = ϖτ(Fu(m), Fv(m), m.constants.ρ₀)
ϖb(m::Model) = ϖb(Fb(m), mixing_depth(m))
dϵ(m::Model, d) = dϵ(ϖτ(m), ϖb(m), d)

unstable(model) = Fb(m) > 0

w_scale_stable(Cκ, Cτb, ϖτ, ϖb, d) = Cκ * ϖτ / (1 + Cτb * d * ϖb^3 / ϖτ^3 )

function w_scale_unstable(Cdϕ, Cτ⁺, Cb⁺, Cτ, Cb, ϖτ, nϕ, dϵ)
    if dϵ < Cdϕ
        return ϖτ * (Cτ⁺ + Cb⁺ * dϵ)^nϕ
    else
        return ϖτ * (Cτ + Cb * dϵ)^(1/3)
    end
end

const nU = 1/4
const nT = 1/2

function w_scale_U(m, i)
    if unstable(m)
        return w_scale_unstable(m.parameters.Cd_U, m.parameters.Cτ⁺_U, m.parameters.Cb⁺_U,
                                m.parameters.Cτ_U, m.parameters.Cb_U,
                                ϖτ(m), nU, dϵ(m, face_depth(m, i))
                                )
    else
        return w_scale_stable(m.parameters.Cκ, m.parameters.Cτb, ϖτ(m), ϖb(m), face_depth(m, i))
    end
end

function w_scale_T(m, i)
    if unstable(m)
        return w_scale_unstable(m.parameters.Cd_T, m.parameters.Cτ⁺_T, m.parameters.Cb⁺_T,
                                m.parameters.Cτ_T, m.parameters.Cb_T,
                                ϖτ(m), nT, dϵ(m, face_depth(m, i))
                                )
    else
        return w_scale_stable(m.parameters.Cκ, m.parameters.Cτb, ϖτ(m), ϖb(m), face_depth(m, i))
    end
end

const w_scale_V = w_scale_U
const w_scale_S = w_scale_T

## ** The K-Profile-Parameterization! **

K_KPP(h, w_scale, d, shape=default_shape_K) = max(0, h * w_scale * shape(d))

# K_{U,V,T,S} is calculated at face points
K_U(m, i) = K_KPP(mixing_depth(m), w_scale_U(m, i), face_depth(m, i)) + m.parameters.K0_U
K_T(m, i) = K_KPP(mixing_depth(m), w_scale_T(m, i), face_depth(m, i)) + m.parameters.K0_T

const K_V = K_U
const K_S = K_T

#
# Interior equations
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
nonlocal_flux(flux, d, shape=default_shape_N) = flux*shape(d) # not minus sign due to flux convention

const N = nonlocal_flux
const d = face_depth

∂NT∂z(m, i) = ( N(FT(m), d(m, i+1)) - N(FT(m), d(m, i)) ) / Δf(m.grid, i)
∂NS∂z(m, i) = ( N(FS(m), d(m, i+1)) - N(FS(m), d(m, i)) ) / Δf(m.grid, i)

# ∇K∇c for c::CellField
const BC = BoundaryCondition

K∂z(K, c, i) = K*∂z(c, i)
∇K∇c(Kᵢ₊₁, Kᵢ, c, i)              = ( K∂z(Kᵢ₊₁, c, i+1) -    K∂z(Kᵢ, c, i)      ) /    Δf(c, i)
∇K∇c_top(Kᵢ, c, top_flux)         = (     -top_flux     - K∂z(Kᵢ, c, length(c)) ) / Δf(c, length(c))
∇K∇c_bottom(Kᵢ₊₁, c, bottom_flux) = (  K∂z(Kᵢ₊₁, c, 2)  +     bottom_flux       ) /    Δf(c, 1)

## Top and bottom flux estimates for constant (Dirichlet) boundary conditions
bottom_flux(K, c, c_bndry, Δf) = -2K*( bottom(c) - c_bndry ) / bottom(Δf) # -K*∂c/∂z at the bottom
top_flux(K, c, c_bndry, Δf)    = -2K*(  c_bndry  -  top(c) ) /   top(Δf)  # -K*∂c/∂z at the top

∇K∇c_top(Kᵢ, c, bc::BC{<:Flux}, model)      = ∇K∇c_top(Kᵢ, c, get_bc(bc, model))
∇K∇c_bottom(Kᵢ₊₁, c, bc::BC{<:Flux}, model) = ∇K∇c_bottom(Kᵢ₊₁, c, get_bc(bc, model))

function calc_rhs_explicit!(rhs, model)
    U, V, T, S = model.solution

    for i in interior(U)
        @inbounds begin
            rhs.U[i] =  f*V[i] + ∇K∇c(K_U(m, i+1), K_U(m, i), U, i)
            rhs.V[i] = -f*U[i] + ∇K∇c(K_V(m, i+1), K_V(m, i), V, i)
            rhs.T[i] =           ∇K∇c(K_T(m, i+1), K_T(m, i), T, i) - ∂NT∂z(m, i)
            rhs.S[i] =           ∇K∇c(K_S(m, i+1), K_S(m, i), S, i) - ∂NS∂z(m, i)
        end
    end

    # Bottom
    rhs.U[1] =  f*bottom(V) * ∇K∇c_bottom(K_U(m, 2), U, model.bcs.U, model)
    rhs.V[1] = -f*bottom(U) * ∇K∇c_bottom(K_V(m, 2), V, model.bcs.V, model)
    rhs.T[1] =                ∇K∇c_bottom(K_T(m, 2), T, model.bcs.T, model) - ∂NT∂z(m, 1)
    rhs.S[1] =                ∇K∇c_bottom(K_S(m, 2), S, model.bcs.S, model) - ∂NS∂z(m, 1)

    # Top
    rhs.U[end] =  f*top(V) * ∇K∇c_top(K_U(m, m.grid.N), U, model.bcs.U, model)
    rhs.V[end] = -f*top(U) * ∇K∇c_top(K_V(m, m.grid.N), V, model.bcs.V, model)
    rhs.T[end] =             ∇K∇c_top(K_T(m, m.grid.N), T, model.bcs.T, model) - ∂NT∂z(m, m.grid.N)
    rhs.S[end] =             ∇K∇c_top(K_S(m, m.grid.N), S, model.bcs.S, model) - ∂NS∂z(m, m.grid.N)

    return nothing
end

end # module
