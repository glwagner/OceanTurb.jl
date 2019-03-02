# Boundary conditions for OceanTurb.jl

struct FieldBoundaryConditions <: FieldVector{2,BoundaryCondition}
  top::BoundaryCondition{:top}
  bottom::BoundaryCondition{:bottom}
end

struct FluxBC{side} <: BoundaryCondition{side}
  flux::Function
end

struct ConstBC{side} <: BoundaryCondition{side}
  value::Function
end
