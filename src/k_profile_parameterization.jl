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
  cP::T
  f::T
end

function Constants(T=Float64; α=2e-4, β=8e-5, ρ₀=1033, cP=3993, f=0) 
  Constants{T}(α, β, ρ₀, cP, f)
end

mutable struct Model{TS,G,T,TP,TC} <: AbstractModel{TS,G,T}
  @add_standard_model_fields
  parameters::Parameters{TP}
  constants::Constants{TC}
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

  return Model(timestepper, grid, solution, equation, bcs, Clock(), parameters, constants)
end

#
# Equation specification
#

@inline α(model) = model.constants.α
@inline β(model) = model.constants.β
@inline T₀(model) = model.constants.T₀
@inline S₀(model) = model.constants.S₀

# Aliases for surface fluxes
@inline FU(model) = model.bcs.U.top.flux(model) 
@inline FV(model) = model.bcs.V.top.flux(model)
@inline FT(model) = model.bcs.T.top.flux(model)
@inline FS(model) = model.bcs.S.top.flux(model)

@inline buoyancy_flux(m) = α(m)*FT(m) + β(m)*FS(m)

## Nonlocal flux

@inline shape_N(d) = 1-d
@inline shape_K(d) = d*(1-d)
@inline nonlocal_flux(flux, d, shape_N=shape_N) = flux*shape_N(d)

@inline ∂NU∂z(model, i) = nonlocal_flux(FU(model), face_depth(model, i))
@inline ∂NV∂z(model, i) = nonlocal_flux(FV(model), face_depth(model, i))
@inline ∂NS∂z(model, i) = nonlocal_flux(FS(model), face_depth(model, i))
@inline ∂NT∂z(model, i) = nonlocal_flux(FT(model), face_depth(model, i))

## Diffusivity

# Mixing depth "h"
@inline mixing_depth(model) = model.grid.Lz # to be changed

# Vertical velocity scale, calculated at face points (insert ref to notes...)
@inline w_scale(Cτ, Cb, Cε, FU, FV, Fb, h, d) = (Cτ*sqrt(FU^2 + FV^2)^3 + Cb*h*min(d*abs(Fb), Cε*Fb))^(1/3)

# Constants can depend on whether field in question is momentum or tracer
@inline w_scale_U(p, m, i) = w_scale(p.Cτ_U, p.Cb_U, p.Cε, FU(m), FV(m),
                                     buoyancy_flux(m), mixing_depth(m), face_depth(m, i))

@inline w_scale_T(p, m, i) = w_scale(p.Cτ_T, p.Cb_T, p.Cε, FU(m), FV(m),
                                     buoyancy_flux(m), mixing_depth(m), face_depth(m, i))

@inline w_scale_U(model, i) = w_scale_U(model.parameters, model, i)
@inline w_scale_T(model, i) = w_scale_T(model.parameters, model, i)

const w_scale_V = w_scale_U
const w_scale_S = w_scale_T

# Non-dimensional depth
@inline face_depth(model, i) = -model.grid.zf[i] / mixing_depth(model)
@inline cell_depth(model, i) = -model.grid.zc[i] / mixing_depth(model)

## ** The K-Profile-Parameterization! **

@inline K_KPP(h, w_scale, d, shape_K=shape_K) = max(0, h*w_scale*shape_K(d))

# K_{U,V,T,S} is calculated at face points
@inline K_U(model, i) = K_KPP(mixing_depth(model), w_scale_U(model, i), face_depth(model, i)) + model.parameters.K0_U
@inline K_V(model, i) = K_KPP(mixing_depth(model), w_scale_V(model, i), face_depth(model, i)) + model.parameters.K0_U
@inline K_T(model, i) = K_KPP(mixing_depth(model), w_scale_T(model, i), face_depth(model, i)) + model.parameters.K0_T
@inline K_S(model, i) = K_KPP(mixing_depth(model), w_scale_S(model, i), face_depth(model, i)) + model.parameters.K0_S

# ∇K∇c for c::CellField
@inline K∂z(K, c, i) = K*∂z(c, i)
@inline ∇K∇c(Kᵢ₊₁, Kᵢ, c::CellField, i)    = ( K∂z(Kᵢ₊₁, c, i+1) -    K∂z(Kᵢ, c, i)  ) /    dzf(c, i)
@inline ∇K∇c_top(K, c::CellField, flux)    = (     -flux      - K∂z(K, c, length(c)) ) / dzf(c, length(c))
@inline ∇K∇c_bottom(K, c::CellField, flux) = (  K∂z(K, c, 2)  +        flux          ) /    dzf(c, 1)


#
# Interior equations
#

@inline ∂U∂t(m, U, V, f, i) =  f*V.data[i] + ∇K∇c(K_U(m, i+1), K_U(m, i), U, i) + ∂NU∂z(m, i)
@inline ∂V∂t(m, U, V, f, i) = -f*U.data[i] + ∇K∇c(K_V(m, i+1), K_V(m, i), U, i) + ∂NV∂z(m, i)
@inline ∂T∂t(m, T, i)       =                ∇K∇c(K_T(m, i+1), K_T(m, i), T, i) + ∂NT∂z(m, i)
@inline ∂S∂t(m, S, i)       =                ∇K∇c(K_S(m, i+1), K_S(m, i), S, i) + ∂NS∂z(m, i)

@inline ∂U∂t(model, i) = ∂U∂t(model, model.solution.U, model.solution.V, model.constants.f, i)
@inline ∂V∂t(model, i) = ∂V∂t(model, model.solution.U, model.solution.V, model.constants.f, i)
@inline ∂T∂t(model, i) = ∂T∂t(model, model.solution.T, i)
@inline ∂S∂t(model, i) = ∂S∂t(model, model.solution.S, i)

#
# Boundary Conditions
#

## Top and bottom flux estimates for constant (Dirichlet) boundary conditions
@inline bottom_flux(K, c, c_bndry, dzf) = -2*K*( bottom(c) - c_bndry ) / bottom(dzf) # -K*∂c/∂z at the bottom
@inline top_flux(K, c, c_bndry, dzf)    = -2*K*(  c_bndry  -  top(c) ) /   top(dzf)  # -K*∂c/∂z at the top

## Flux BCs --- these are incorrect and must be fixed.
@inline ∂U∂t(model, bc::FluxBC{Top}) = 0.0 #∇K∇c_top(K_U(model, model.grid.nz), model.solution.U, bc.flux(model))
@inline ∂V∂t(model, bc::FluxBC{Top}) = 0.0 #∇K∇c_top(K_V(model, model.grid.nz), model.solution.V, bc.flux(model))
@inline ∂T∂t(model, bc::FluxBC{Top}) = 0.0 #∇K∇c_top(K_T(model, model.grid.nz), model.solution.T, bc.flux(model))
@inline ∂S∂t(model, bc::FluxBC{Top}) = 0.0 #∇K∇c_top(K_S(model, model.grid.nz), model.solution.S, bc.flux(model))

@inline ∂U∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_U(model, 2), model.solution.U, bc.flux(model))
@inline ∂V∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_V(model, 2), model.solution.V, bc.flux(model))
@inline ∂T∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_T(model, 2), model.solution.T, bc.flux(model))
@inline ∂S∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_S(model, 2), model.solution.S, bc.flux(model))

## Value BCs
# ** note these are incorrect because neither Coriolis terms nor non-local flux are included **
function ∂U∂t(model, bc::ValueBC{Bottom}) 
  flux = bottom_flux(K_U(model, 1), model.solution.U, bc.value(model), model.grid.dzf)
  return ∇K∇c_bottom(K_U(model, 2), model.solution.U, flux)
end

function ∂V∂t(model, bc::ValueBC{Bottom}) 
  flux = bottom_flux(K_V(model, 1), model.solution.V, bc.value(model), model.grid.dzf)
  return ∇K∇c_bottom(K_V(model, 2), model.solution.V, flux)
end

function ∂T∂t(model, bc::ValueBC{Bottom}) 
  flux = bottom_flux(K_T(model, 1), model.solution.T, bc.value(model), model.grid.dzf)
  return ∇K∇c_bottom(K_T(model, 2), model.solution.T, flux)
end

function ∂S∂t(model, bc::ValueBC{Bottom}) 
  flux = bottom_flux(K_S(model, 1), model.solution.S, bc.value(model), model.grid.dzf)
  return ∇K∇c_bottom(K_S(model, 2), model.solution.S, flux)
end

function ∂U∂t(model, bc::ValueBC{Top}) 
  flux = top_flux(K_U(model, model.grid.nz+1), model.solution.U, bc.value(model), model.grid.dzf)
  return ∇K∇c_top(K_U(model, model.grid.nz  ), model.solution.U, flux)
end

function ∂V∂t(model, bc::ValueBC{Top}) 
  flux = top_flux(K_V(model, model.grid.nz+1), model.solution.V, bc.value(model), model.grid.dzf)
  return ∇K∇c_top(K_V(model, model.grid.nz  ), model.solution.V, flux)
end

function ∂T∂t(model, bc::ValueBC{Top}) 
  flux = top_flux(K_T(model, model.grid.nz+1), model.solution.T, bc.value(model), model.grid.dzf)
  return ∇K∇c_top(K_T(model, model.grid.nz  ), model.solution.T, flux)
end

function ∂S∂t(model, bc::ValueBC{Top}) 
  flux = top_flux(K_S(model, model.grid.nz+1), model.solution.S, bc.value(model), model.grid.dzf)
  return ∇K∇c_top(K_S(model, model.grid.nz  ), model.solution.S, flux)
end

end # module
