# Boundary conditions for OceanTurb.jl

@inline zero_function(args...) = 0

struct FieldBoundaryConditions <: FieldVector{2,BoundaryCondition}
  top::BoundaryCondition{:top}
  bottom::BoundaryCondition{:bottom}
end

function FieldBoundaryConditions(;
       top = FluxBC{:top}(zero_function), 
    bottom = FluxBC{:bottom}(zero_function)
  )
  FieldBoundaryConditions(top, bottom)
end

struct FluxBC{side} <: BoundaryCondition{side}
  flux::Function
end

struct ConstBC{side} <: BoundaryCondition{side}
  value::Function
end
