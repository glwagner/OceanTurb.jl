module Diffusion

using 
  Reexport

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

DκDc(κ::Number, c, i) = κ * ∂²z(c, i)
DκDc(κ::Field, c, i) = (κ.data[i+1]*∂z(c, i+1) - κ.data[i]*∂z(c, i)) / dzf(c, i)

function dcdt!(rhs, model)
  for i in eachindex(rhs)
    @inbounds rhs.data[i] = DκDc(model.parameters.κ, model.solution, i)
  end
  nothing
end 

mutable struct Model{G,S,TS,T} <: AbstractModel{G,S,TS}  
  grid::G
  solution::S
  timestepper::TS
  parameters::Parameters
  clock::Clock{T}
end

function Model(;
       nz = 100, 
       Lz = 1.0,
        κ = 0.1,
  stepper = :ForwardEuler
)

  grid = UniformGrid(nz, Lz)
  parameters = Parameters(κ)
  solution = Solution(CellField(grid))
  timestepper = Timestepper(:ForwardEuler, solution)

  Model(grid, solution, timestepper, parameters, Clock())
end

function set_kappa!(model, κ)
  model.parameters = Parameters(κ)
  return nothing
end

# Forward Euler timestepping
function step!(model, dt)
  dcdt!(model.timestepper.rhs, model)
  @. model.solution.c.data += dt*model.timestepper.rhs.data
  return nothing
end

function step!(model, dt, nt)
  for step = 1:nt
    step!(model, dt)
  end
  return nothing
end

end # module
