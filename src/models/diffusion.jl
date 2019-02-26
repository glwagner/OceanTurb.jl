module Diffusion

using OceanTurb

export
  Model,
  set!

struct DiffusionParameters{T} <: Parameters
  κ::T
end

function diffusion_operator!(rhs, model)
  for i in eachindex(rhs)
    @inbounds rhs[i] = DκDc(i, model.parameters.κ, model.solution)
  end
  nothing
end 

DκDc(κ::Number, c, i) = κ * ∂²z(c, i)
DκDc(κ::Field, c, i) = (κ.data[i+1]*∂z(c, i+1) - κ.data[i]*∂z(c, i)) / dzf(c, i)

function Model(;
  nz = 100, 
  Lz = 1.0,
   κ = 0.1
  )

  grid = UniformGrid(nz, Lz)
  parameters = DiffusionParameters(κ)
  equation = StandardEquation(nothing, diffusion_operator!)
  constants = nothing
  solution = CellField(grid)

  OceanTurb.Model(grid, constants, parameters, equation, solution)
end

function set!(model, c)
  for i in eachindex(c)
    @inbounds model.solution.data[i] = c[i]
  end
  return nothing
end

set!(model; c=model.solution) = set!(model, c)

end # module
