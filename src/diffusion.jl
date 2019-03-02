module Diffusion

using 
  Reexport,
  StaticArrays

@reexport using OceanTurb

export
  Model,
  step!

struct Parameters{T} <: AbstractParameters
  κ::T
end

struct Solution <: AbstractSolution{1,CellField}
  c::CellField
end

abstract type BoundaryCondition{side} end

struct BoundaryConditions <: FieldVector{2,BoundaryCondition}
  top::BoundaryCondition{:top}
  bottom::BoundaryCondition{:bottom}
end

@inline zero_function(args...) = 0

struct FluxBC{side} <: BoundaryCondition{side}
  flux::Function
end

struct ConstBC{side} <: BoundaryCondition{side}
  value::Function
end

mutable struct Model{G,S,TS,T,B} <: AbstractModel{G,S,TS}  
  grid::G
  solution::S
  timestepper::TS
  parameters::Parameters
  bcs::B
  clock::Clock{T}
end

function Model(;
       nz = 100, 
       Lz = 1.0,
        κ = 0.1,
  stepper = :ForwardEuler,
  bcs = BoundaryConditions(FluxBC{:top}(zero_function), FluxBC{:bottom}(zero_function))
)

  grid = UniformGrid(nz, Lz)
  parameters = Parameters(κ)
  solution = Solution(CellField(grid))
  timestepper = Timestepper(:ForwardEuler, solution)

  Model(grid, solution, timestepper, parameters, bcs, Clock())
end

function set_kappa!(model, κ)
  model.parameters = Parameters(κ)
  return nothing
end

# Operators for the interior and boundary
κ∂z(κ::Number, c, i) = κ*∂z(c, i)
κ∂z(κ::Field, c, i) = κ.data[i]*∂z(c, i)

@inline ∇κ∇c(κ::Number, c, i)   = κ * ∂²z(c, i)
@inline ∇κ∇c(κ, c, i)           = ( κ∂z(κ, c, i+1) -    κ∂z(κ, c, i)      ) /    dzf(c, i)
@inline ∇κ∇c_top(κ, c, flux)    = (     -flux      - κ∂z(κ, c, length(c)) ) / dzf(c, length(c))
@inline ∇κ∇c_bottom(κ, c, flux) = (  κ∂z(κ, c, 2)  +        flux          ) /    dzf(c, 1)

# Interior
@inline ∂ₜc(model, i) = ∇κ∇c(model.parameters.κ, model.solution.c, i)

# Boundary conditions
∂ₜc(model, bc::BoundaryCondition) = throw("Boundary condition $bc not implemented.")
∂ₜc(model, bc::FluxBC{:top}) = ∇κ∇c_top(model.parameters.κ, model.solution.c, bc.flux(model))
∂ₜc(model, bc::FluxBC{:bottom}) = ∇κ∇c_bottom(model.parameters.κ, model.solution.c, bc.flux(model))

# Top and bottom flux estimates for constant (Dirichlet) boundary conditions
@inline bottom_flux(κ, c, c_boundary, dzf) = -2*bottom(κ)*( c_boundary - bottom(c)  ) / bottom(dzf) # -κ*∂c/∂z at the bottom
@inline top_flux(κ, c, c_boundary, dzf)    = -2*  top(κ) *(   top(c)   - c_boundary ) /   top(dzf)  # -κ*∂c/∂z at the top

function ∂ₜc(model, bc::ConstBC{:bottom}) 
  flux = bottom_flux(model.parameters.κ, model.solution.c, bc.value(model), model.grid.dzf)
  return ∇κ∇c_bottom(model.parameters.κ, model.solution.c, flux)
end

function ∂ₜc(model, bc::ConstBC{:top}) 
  flux = top_flux(model.parameters.κ, model.solution.c, bc.value(model), model.grid.dzf)
  return ∇κ∇c_top(model.parameters.κ, model.solution.c, flux)
end

# Forward Euler timestepping
function step!(model, Δt)

  # Interior step
  for i in interior(model.solution.c)
    @inbounds model.timestepper.rhs.c.data[i] = ∂ₜc(model, i)
  end

  # Boundary conditions
  model.timestepper.rhs.c.data[end] = ∂ₜc(model, model.bcs.top)
  model.timestepper.rhs.c.data[1]   = ∂ₜc(model, model.bcs.bottom)

  # Time step
  @. model.solution.c.data += Δt*model.timestepper.rhs.c.data

  return nothing
end

#= 
Sketch of RK2 time-stepping...
interior:
RHS[i] = ∂ₜc(model, i)

boundary:
RHS[1]   = ∂ₜc(model, model.bcs.bottom)
RHS[end] = ∂ₜc(model, model.bcs.top)

etc.
=#

function step!(model, Δt, nt)
  for step = 1:nt
    step!(model, Δt)
  end
  return nothing
end

end # module
