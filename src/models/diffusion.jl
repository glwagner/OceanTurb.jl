module Diffusion

using OceanTurb

struct DiffusionParameters{T} <: Parameters
  κ::T
end

function diffusion_operator!(rhs, model)
  for i in eachindex(rhs)
    @inbounds rhs[i] = DκDc(i, model.parameters.κ, model.solution)
  end
  nothing
end 

DκDc(κ::Number, c, i) = κ * (∂z(i+1, c) - ∂z(i, c)) / dzf(c)
DκDc(κ::Field, c, i) = (κ.data[i+1]*∂z(c, i+1) - κ.data[i]*∂z(c, i)) / dzf(c, i)

function Model(;
  nz = 100, 
  Lz = 1.0
   κ = 0.1,
  )

  grid = UniformGrid(nz, Lz)
  parameters = DiffusionParameters(κ)
  equation = StandardEquation(NullOperator, diffusion_operator!)
  constants = nothing
  solution = CellField(grid)

  Model(grid, constants, parameters, equation, solution)
end

end # module
