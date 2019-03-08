module KPP

using
  Reexport,
  StaticArrays

@reexport using OceanTurb

export
  Parameters,
  Constants,
  Model

const nsol = 4

@specify_solution CellField U V T S

struct Parameters{T} <: AbstractParameters
  Cε::T    # Surface layer fraction
  Cτ_U::T  # Effect of surface momentum flux on momentum diffusivity
  Cb_U::T  # Effect of surface buoyancy flux on momentum diffusivity
  Cτ_T::T  # Effect of surface momentum flux on momentum diffusivity
  Cb_T::T  # Effect of surface buoyancy flux on momentum diffusivity
  CRi::T   # Critical Richardson number
  CKE::T   # Unresolved turbulence parameter
  K0_U::T
  K0_T::T
  K0_S::T
end

function Parameters(T=Float64; Cε=0.1, Cτ_U=1, Cb_U=1, Cτ_T=1, Cb_T=1, CRi=1, CKE=1,
                    K0_U=1e-5, K0_T=1e-5, K0_S=1e-5)
  Parameters{T}(Cε, Cτ_U, Cb_U, Cτ_T, Cb_T, CRi, CKE, K0_U, K0_T, K0_S)
end

struct Constants{T}
  α::T
  β::T
  ρ₀::T
  g::T
  cP::T
  f::T
end

function Constants(T=Float64; α=2e-4, β=8e-5, ρ₀=1033, cP=3993, f=0, g=9.81)
  Constants{T}(α, β, ρ₀, cP, f)
end

mutable struct Model{TS, G, T, TP, TC, A} <: AbstractModel{TS, G, T}
  @add_standard_model_fields
  parameters::Parameters{TP}
  constants::Constants{TC}
  δRi::FaceField{A, G}
end

function Model(;
              nz = 100,
              Lz = 1.0,
              K0 = 1e-5, K0_U = K0, K0_T = K0, K0_S = K0,
         stepper = :ForwardEuler,
             bcs = BoundaryConditions((ZeroFlux() for i=1:nsol)...),
       constants = Constants(),
      parameters = Parameters(K0_U=K0_U, K0_T=K0_T, K0_S=K0_S)
)

  grid = UniformGrid(nz, Lz)
  solution = Solution((CellField(grid) for i=1:nsol)...)
  equation = Equation(∂U∂t, ∂V∂t, ∂T∂t, ∂S∂t)
  timestepper = Timestepper(:ForwardEuler, solution)
  δRi = FaceField(grid)

  return Model(timestepper, grid, solution, equation, bcs, Clock(), parameters, constants, δRi)
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
FU(m) = m.bcs.U.top.flux(m)
FV(m) = m.bcs.V.top.flux(m)
FT(m) = m.bcs.T.top.flux(m)
FS(m) = m.bcs.S.top.flux(m)
Fb(m) = g(m)/ρ₀(m) * (α(m)*FT(m) - β(m)*FS(m))

const buoyancy_flux = Fb
Bz(α, β, T, S, i) = α*∂z(T, i) - β*∂z(S, i)

# Shape function. 'd' is a non-dimensional depth coordinate.
default_shape_N(d) = d*(1-d)^2
default_shape_K(d) = d*(1-d)^2

"""
    nonlocal_flux(flux, d, shape)

Returns the nonlocal flux, N = -flux*shape(d),
where `flux` is the flux of some quantity out of the surface,
`shape` is a shape function, and `d` is a non-dimensional depth coordinate
that increases from 0 at the surface to 1 at the bottom of the mixing layer.

Because flux is defined as pointing in the positive direction,
a positive surface flux implies negative surface flux divergence,
which implies a reduction to the quantity in question.
For example, positive heat flux out of the surface implies cooling.
"""
nonlocal_flux(flux, d, shape=default_shape_N) = flux*shape(d) # not minus sign due to flux convention

∂NU∂z(m, i) = nonlocal_flux(FU(m), face_depth(m, i))
∂NV∂z(m, i) = nonlocal_flux(FV(m), face_depth(m, i))
∂NS∂z(m, i) = nonlocal_flux(FS(m), face_depth(m, i))
∂NT∂z(m, i) = nonlocal_flux(FT(m), face_depth(m, i))

## Diffusivity

# Mixing depth "h"

"Returns the surface_layer_average for mixing depth h=- zf[i]."
function surface_layer_average(Cε, c, i)
  iε = 1 + Cε * (length(c) + 1 - i) # (fractional) "index" of the surface layer

  # Contribution of fractional cell to total integral
  surface_layer_integral = frac * dzf(c, face-1) * c.data[face-1]

  # Add cells above face, if there are any.
  for j = face:length(c)
    @inbounds surface_layer_integral += dzf(c, j) * c.data[j]
  end

  h = -c.grid.zf[i] # depth
  return surface_layer_integral / (Cε*h)
end

"""
Return Δc(hᵢ), the difference between the surface-layer average of c and its value at depth hᵢ, where
i is a face index.
"""
Delta(Cε, c, i) = surface_layer_average(Cε, c, i) - avz(c, i)

"Returns the parameterization for unresolved KE at face point i."
function unresolved_KE(CKE, Fb, α, β, T, S, i)
  h = -T.grid.zf[i]
  return CKE * h^(4/3) * sqrt(Bz(α, β, T, S, i)) * max(0, Fb^(1/3))
end

"""
    Richardson(model, i)

Returns the Richardson number of `model`, including an estimate
for unresolved kinetic energy, at face point i.
"""
function Richardson(CKE, Cε, Fb, α, β, g, ρ₀, U, V, T, S, i)
  h = -T.grid.zf[i]
   ΔB = g/ρ₀ * (α*Delta(Cε, T, i) - β*Delta(Cε, S, i))
  ΔU² = Delta(Cε, U, i)^2 + Delta(Cε, V, i)^2
  return h*ΔB / (ΔU² + unresolved_KE(CKE, Fb, α, β, T, S, i))
end

Richardson(m, i) = Richardson(m.parameters.CKE, m.parameters.Cε, buoyancy_flux(m),
                              α(m), β(m), g(m), ρ₀(m), m.solution.U, m.solution.V,
                              m.solution.T, m.solution.S, i)

"""
    mixing_depth(model)

Calculate 'h', the mixing depth, for the current state of `model`.
This function calculates δRi.
"""
function mixing_depth(m)
  # Calculate deviation of Ri from critical value at each face point
  for i = eachindex(m.δRi)
    m.δRi.data[i] = m.parameters.CRi - Richardson(m, i)
  end
  # Linearly interpolate to find where δRi = 0.
  ih₂ = searchsortedfirst(m.δRi.data, 0)
  ih₁ = max(1, ih₂ - 1) # don't let ih₁ < 1.
  Δz = m.grid.zf[ih₂] - m.grid.zf[ih₁]
  ΔδRi = m.δRi.data[ih₂] - m.δRi.data[ih₁] # Delta of Ri deviation
  z★ = m.grid.zf[ih₂] - m.δRi.data[ih₂] * Δz / ΔδRi
  return -z★
end

# Vertical velocity scale, calculated at face points (insert ref to notes...)
w_scale(Cτ, Cb, Cε, FU, FV, Fb, h, d) = (Cτ*sqrt(FU^2 + FV^2)^3 + Cb*h*min(d*abs(Fb), Cε*Fb))^(1/3)

# Parameters can depend on whether field in question is momentum or tracer
w_scale_U(p, m, i) = w_scale(p.Cτ_U, p.Cb_U, p.Cε, FU(m), FV(m),
                                     buoyancy_flux(m), mixing_depth(m), face_depth(m, i))

w_scale_T(p, m, i) = w_scale(p.Cτ_T, p.Cb_T, p.Cε, FU(m), FV(m),
                                     buoyancy_flux(m), mixing_depth(m), face_depth(m, i))

w_scale_U(m, i) = w_scale_U(m.parameters, m, i)
w_scale_T(m, i) = w_scale_T(m.parameters, m, i)

const w_scale_V = w_scale_U
const w_scale_S = w_scale_T

# Non-dimensional depth
face_depth(m, i) = -m.grid.zf[i] / mixing_depth(m)
cell_depth(m, i) = -m.grid.zc[i] / mixing_depth(m)

## ** The K-Profile-Parameterization! **

K_KPP(h, w_scale, d, shape=default_shape_K) = max(0, h*w_scale*shape(d))

# K_{U,V,T,S} is calculated at face points
K_U(m, i) = K_KPP(mixing_depth(m), w_scale_U(m, i), face_depth(m, i)) + m.parameters.K0_U
K_V(m, i) = K_KPP(mixing_depth(m), w_scale_V(m, i), face_depth(m, i)) + m.parameters.K0_U
K_T(m, i) = K_KPP(mixing_depth(m), w_scale_T(m, i), face_depth(m, i)) + m.parameters.K0_T
K_S(m, i) = K_KPP(mixing_depth(m), w_scale_S(m, i), face_depth(m, i)) + m.parameters.K0_S

# ∇K∇c for c::CellField
const CF = CellField
const FF = FaceField

K∂z(K, c, i) = K*∂z(c, i)
∇K∇c(Kᵢ₊₁, Kᵢ, c::CF, i)              = ( K∂z(Kᵢ₊₁, c, i+1) -    K∂z(Kᵢ, c, i)      ) /    dzf(c, i)
∇K∇c_top(Kᵢ, c::CF, top_flux)         = (     -top_flux     - K∂z(Kᵢ, c, length(c)) ) / dzf(c, length(c))
∇K∇c_bottom(Kᵢ₊₁, c::CF, bottom_flux) = (  K∂z(Kᵢ₊₁, c, 2)  +     bottom_flux       ) /    dzf(c, 1)

#
# Interior equations
#

∂U∂t(m, U, V, f, i) =  f*V.data[i] + ∇K∇c(K_U(m, i+1), K_U(m, i), U, i) - ∂NU∂z(m, i)
∂V∂t(m, U, V, f, i) = -f*U.data[i] + ∇K∇c(K_V(m, i+1), K_V(m, i), U, i) - ∂NV∂z(m, i)
∂T∂t(m, T, i)       =                ∇K∇c(K_T(m, i+1), K_T(m, i), T, i) - ∂NT∂z(m, i)
∂S∂t(m, S, i)       =                ∇K∇c(K_S(m, i+1), K_S(m, i), S, i) - ∂NS∂z(m, i)

∂U∂t(m, i) = ∂U∂t(m, m.solution.U, m.solution.V, m.constants.f, i)
∂V∂t(m, i) = ∂V∂t(m, m.solution.U, m.solution.V, m.constants.f, i)
∂T∂t(m, i) = ∂T∂t(m, m.solution.T, i)
∂S∂t(m, i) = ∂S∂t(m, m.solution.S, i)

#
# Boundary Conditions specification
#

## Top and bottom flux estimates for constant (Dirichlet) boundary conditions
bottom_flux(K, c, c_bndry, dzf) = -2*K*( bottom(c) - c_bndry ) / bottom(dzf) # -K*∂c/∂z at the bottom
top_flux(K, c, c_bndry, dzf)    = -2*K*(  c_bndry  -  top(c) ) /   top(dzf)  # -K*∂c/∂z at the top

## Flux BCs --- omit diffusive flux at the top for KPP (?)
∂U∂t(m, U, V, f, ::FluxBC{Top}) =  f*top(V) - ∂NU∂z(m, length(U))
∂V∂t(m, U, V, f, ::FluxBC{Top}) = -f*top(U) - ∂NV∂z(m, length(V))
∂T∂t(m, T, ::FluxBC{Top})       =             ∂NT∂z(m, length(T))
∂S∂t(m, S, ::FluxBC{Top})       =             ∂NS∂z(m, length(S))

∂U∂t(m, U, V, f, bc::FluxBC{Bottom}) =  f*bottom(V) + ∇K∇c_bottom(K_U(m, 2), U, bc.flux(m)) - ∂NU∂z(m, 1)
∂V∂t(m, U, V, f, bc::FluxBC{Bottom}) = -f*bottom(U) + ∇K∇c_bottom(K_V(m, 2), V, bc.flux(m)) - ∂NV∂z(m, 1)
∂T∂t(m, T, bc::FluxBC{Bottom})       =                ∇K∇c_bottom(K_T(m, 2), T, bc.flux(m)) - ∂NT∂z(m, 1)
∂S∂t(m, S, bc::FluxBC{Bottom})       =                ∇K∇c_bottom(K_S(m, 2), S, bc.flux(m)) - ∂NS∂z(m, 1)

∂U∂t(m, bc::FluxBC{Bottom}) = ∂U∂t(m, m.solution.U, m.solution.V, m.constants.f, bc)
∂V∂t(m, bc::FluxBC{Bottom}) = ∂V∂t(m, m.solution.U, m.solution.V, m.constants.f, bc)
∂T∂t(m, bc::FluxBC{Bottom}) = ∂T∂t(m, m.solution.T, bc)
∂S∂t(m, bc::FluxBC{Bottom}) = ∂S∂t(m, m.solution.S, bc)

## Value BCs -- implemented at the bottom only.
∂U∂t(m, bc::ValueBC{Top}) = throw("A top ValueBC cannot be defined in KPP.")
∂V∂t(m, bc::ValueBC{Top}) = throw("A top ValueBC cannot be defined in KPP.")
∂T∂t(m, bc::ValueBC{Top}) = throw("A top ValueBC cannot be defined in KPP.")
∂S∂t(m, bc::ValueBC{Top}) = throw("A top ValueBC cannot be defined in KPP.")

∂U∂t(V, f, ∇K∇U, ∂N∂z, bc::ValueBC{Bottom}) = -f*V + ∇K∇U - ∂N∂z
∂V∂t(U, f, ∇K∇V, ∂N∂z, bc::ValueBC{Bottom}) =  f*U + ∇K∇V - ∂N∂z
∂T∂t(      ∇K∇T, ∂N∂z, bc::ValueBC{Bottom}) = ∇K∇T + ∂N∂z
∂S∂t(      ∇K∇S, ∂N∂z, bc::ValueBC{Bottom}) = ∇K∇S + ∂N∂z

function ∂U∂t(m, bc::ValueBC{Bottom})
  flux = bottom_flux(K_U(m, 1), U, bc.value(m), m.grid.dzf)
  return ∂U∂t(bottom(m.solution.V), m.constants.f, ∇K∇c_bottom(K_U(m, 2), U, flux), ∂NU∂z(m, 1))
end

function ∂V∂t(m, bc::ValueBC{Bottom})
  flux = bottom_flux(K_V(m, 1), U, bc.value(m), m.grid.dzf)
  return ∂V∂t(bottom(m.solution.U), m.constants.f, ∇K∇c_bottom(K_V(m, 2), V, flux), ∂NV∂z(m, 1))
end

function ∂T∂t(m, bc::ValueBC{Bottom})
  flux = bottom_flux(K_T(m, 1), T, bc.value(m), m.grid.dzf)
  return ∂T∂t(∇K∇c_bottom(K_T(m, 2), T, flux), ∂NT∂z(m, 1))
end

function ∂S∂t(m, bc::ValueBC{Bottom})
  flux = bottom_flux(K_S(m, 1), S, bc.value(m), m.grid.dzf)
  return ∂S∂t(∇K∇c_bottom(K_S(m, 2), S, flux), ∂NS∂z(m, 1))
end

end # module
