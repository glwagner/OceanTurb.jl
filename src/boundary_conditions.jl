# Boundary conditions for OceanTurb.jl

abstract type Boundary end
struct Top <: Boundary end
struct Bottom <: Boundary end

abstract type BoundaryCondition{B<:Boundary} end

@inline zero_function(args...) = 0

struct FieldBoundaryConditions <: FieldVector{2,BoundaryCondition}
  bottom::BoundaryCondition{Bottom}
  top::BoundaryCondition{Top}
end

function FieldBoundaryConditions(;
    bottom = FluxBC{Bottom}(zero_function),
       top = FluxBC{Top}(zero_function),
  )
  FieldBoundaryConditions(bottom, top)
end

ZeroFlux() = FieldBoundaryConditions() # the default

struct FluxBC{B<:Boundary} <: BoundaryCondition{B}
  flux::Function
end

struct ValueBC{B<:Boundary} <: BoundaryCondition{B}
  value::Function
end
