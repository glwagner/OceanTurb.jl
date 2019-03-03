module KPP

using 
  Reexport,
  StaticArrays

@reexport using OceanTurb

export
  Parameters,
  SurfaceFluxes,
  Constants,
  Model,
  step!

@specify_solution CellField U V T S

struct Parameters <: AbstractParameters
  const_ε::Float64    # Surface layer fraction
  const_U_u::Float64  # Effect of surface momentum flux on momentum diffusivity
  const_U_b::Float64  # Effect of surface buoyancy flux on momentum diffusivity
  const_C_u::Float64  # Effect of surface momentum flux on scalar diffusivity
  const_C_b::Float64  # Effect of surface buoyancy flux on scalar diffusivity
  const_C_Ri::Float64 # Critical Richardson number
  const_C_KE::Float64 # Unresolved turbulence parameter
end

struct SurfaceFluxes{TT,TS,TU}
  FT::TT
  FS::TS
  FU::TU
  FV::TU
end

struct Constants
  α::Float64
  β::Float64
  ρ₀::Float64
  cP::Float64
  K0::Float64
  function Constants(α=2e-4, β=8e-5, ρ₀=1033, cP=3993, K0=1e-5) 
    new(α, β, ρ₀, cP, K0)
  end
end


mutable struct Model{FT,FS,FU,TS,G,T} <: AbstractModel{TS,G,T}
  @add_standard_model_fields
  parameters::Parameters
  surface_fluxes::SurfaceFluxes{FT,FS,FU}
  constants::Constants
end

function Model(;
              nz = 100, 
              Lz = 1.0,
         stepper = :ForwardEuler,
             bcs = BoundaryConditions(ZeroFlux(), ZeroFlux(), ZeroFlux(), ZeroFlux()),
  surface_fluxes = SurfaceFluxes(0, 0, 0, 0),
       constants = Constants(),
      parameters = Parameters(1, 1, 1, 1, 1, 1, 1, 1)
)

  grid = UniformGrid(nz, Lz)
  solution = Solution(CellField(grid))
  equation = Equation(∂c∂t)
  timestepper = Timestepper(:ForwardEuler, solution)

  return Model(timestepper, grid, solution, equation, bcs, Clock(), parameters)
end

#
# Equation specification
#
@inline α(model) = model.constants.α
@inline β(model) = model.constants.β
@inline T₀(model) = model.constants.T₀
@inline S₀(model) = model.constants.S₀
@inline FU(model) = model.surface_fluxes.FU(model)
@inline FV(model) = model.surface_fluxes.FV(model)
@inline FT(model) = model.surface_fluxes.FT(model)
@inline FS(model) = model.surface_fluxes.FS(model)

@inline buoyancy_flux(model) = α(model)*FT(model) + β(model)*FS(model)

## Nonlocal flux

@inline nonlocal_flux(const_N, flux, d) = const_N*flux*(1-d)
@inline nonlocal_flux(flux, d) = const_N*flux*(1-d)

@inline ∂NU∂z(model, i) = nonlocal_flux(model.surface_fluxes.FU(model), face_depth(model, i))
@inline ∂NV∂z(model, i) = nonlocal_flux(model.surface_fluxes.FV(model), face_depth(model, i))
@inline ∂NS∂z(model, i) = nonlocal_flux(model.surface_fluxes.FS(model), face_depth(model, i))
@inline ∂NT∂z(model, i) = nonlocal_flux(model.surface_fluxes.FT(model), face_depth(model, i))

## Diffusivity: depth and vertical velocity scale

mixing_depth(model) = model.grid.Lz # to be changed

@inline face_depth(model, i) = -model.grid.zf[i] / mixing_depth(model)
@inline cell_depth(model, i) = -model.grid.zc[i] / mixing_depth(model)

@inline function w_scale(const_u, const_b, const_ε, abs_momentum_flux, buoyancy_flux, h, d)
  (   const_u * sqrt(abs_momentum_flux)^3 
    + const_b * h * min(d*abs(buoyancy_flux), const_ε*buoyancy_flux) )^(1/3)
end

@inline function w_scale_U(model, i)
  w_scale(model.parameters.const_U_u, model.parameters.const_U_b, model.parameters.const_ε, 
          sqrt(FU(model)^2 + FV(model)^2), buoyancy_flux(model), mixing_depth(model), face_depth(model, i))
end

@inline function w_scale_T(model, i)
  w_scale(model.parameters.const_C_u, model.parameters.const_C_b, model.parameters.const_ε, 
          sqrt(FU(model)^2 + FV(model)^2), buoyancy_flux(model), mixing_depth(model), face_depth(model, i))
end

@inline w_scale_V(model, i) = w_scale_U(model, i)
@inline w_scale_S(model, i) = w_scale_T(model, i)

# ** KPP! **
@inline K_KPP(h, w_scale, d) = max(0, h * w_scale * d * (1 - d)^2)

@inline K_U(model, i) = K_KPP(mixing_depth(model), w_scale_U(model, i), face_depth(model, i)) + model.constants.K0
@inline K_V(model, i) = K_KPP(mixing_depth(model), w_scale_V(model, i), face_depth(model, i)) + model.constants.K0
@inline K_T(model, i) = K_KPP(mixing_depth(model), w_scale_T(model, i), face_depth(model, i)) + model.constants.K0
@inline K_S(model, i) = K_KPP(mixing_depth(model), w_scale_S(model, i), face_depth(model, i)) + model.constants.K0

K∂z(K::Number, c, i) = K*∂z(c, i)
K∂z(K::AbstractField, c, i) = K.data[i]*∂z(c, i)
K∂z(K::Function, c, i) = K(c.grid.zf[i]) * ∂z(c, i) # works for K(z)

@inline ∇K∇c(Kᵢ₊₁, Kᵢ, c, i)    = ( K∂z(Kᵢ₊₁, c, i+1) -    K∂z(Kᵢ, c, i)  ) /    dzf(c, i)
@inline ∇K∇c_top(K, c, flux)    = (     -flux      - K∂z(K, c, length(c)) ) / dzf(c, length(c))
@inline ∇K∇c_bottom(K, c, flux) = (  K∂z(K, c, 2)  +        flux          ) /    dzf(c, 1)


#
# Interior equations
#

@inline ∂U∂t(model, i) = ∇K∇c(K_U(model, i+1), K_U(model, i), model.solution.U, i) + ∂NU∂z(model, i)
@inline ∂V∂t(model, i) = ∇K∇c(K_V(model, i+1), K_V(model, i), model.solution.V, i) + ∂NV∂z(model, i)
@inline ∂T∂t(model, i) = ∇K∇c(K_T(model, i+1), K_T(model, i), model.solution.S, i) + ∂NT∂z(model, i)
@inline ∂S∂t(model, i) = ∇K∇c(K_S(model, i+1), K_S(model, i), model.solution.T, i) + ∂NS∂z(model, i)
 

#
# Boundary Conditions
#

## Top and bottom flux estimates for constant (Dirichlet) boundary conditions
@inline bottom_flux(K, c, c_bndry, dzf) = -2*K*( bottom(c) - c_bndry ) / bottom(dzf) # -K*∂c/∂z at the bottom
@inline top_flux(K, c, c_bndry, dzf)    = -2*K*(  c_bndry  -  top(c) ) /   top(dzf)  # -K*∂c/∂z at the top

## Flux BCs
@inline ∂U∂t(model, bc::FluxBC{Top})    = ∇K∇c_top(K_U(model, model.grid.nz), model.solution.U, bc.flux(model))
@inline ∂V∂t(model, bc::FluxBC{Top})    = ∇K∇c_top(K_V(model, model.grid.nz), model.solution.V, bc.flux(model))
@inline ∂T∂t(model, bc::FluxBC{Top})    = ∇K∇c_top(K_T(model, model.grid.nz), model.solution.T, bc.flux(model))
@inline ∂S∂t(model, bc::FluxBC{Top})    = ∇K∇c_top(K_S(model, model.grid.nz), model.solution.S, bc.flux(model))

@inline ∂U∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_U(model, 2), model.solution.U, bc.flux(model))
@inline ∂V∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_V(model, 2), model.solution.V, bc.flux(model))
@inline ∂T∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_T(model, 2), model.solution.T, bc.flux(model))
@inline ∂S∂t(model, bc::FluxBC{Bottom}) = ∇K∇c_bottom(K_S(model, 2), model.solution.S, bc.flux(model))

## Value BCs
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
