# Boundary conditions for OceanTurb.jl

import Base: convert

abstract type Boundary end
struct Top <: Boundary end
struct Bottom <: Boundary end
struct Unset <: Boundary end

abstract type BoundaryCondition{B<:Boundary} end

@inline zero_function(args...) = 0

mutable struct FieldBoundaryConditions <: FieldVector{2,BoundaryCondition}
  bottom::BoundaryCondition{Bottom}
  top::BoundaryCondition{Top}
end

struct FluxBC{B<:Boundary} <: BoundaryCondition{B}
  flux::Function
end

struct ValueBC{B<:Boundary} <: BoundaryCondition{B}
  value::Function
end

# 
# Sugary goodness
#

convert(::Type{FluxBC{B}}, flux_bc::FluxBC{C}) where {B,C} = FluxBC{B}(flux_bc.flux)
convert(::Type{ValueBC{B}}, value_bc::ValueBC{C}) where {B,C} = ValueBC{B}(value_bc.value)

function FluxBC(boundary, flux::Number)
  @inline flux_function(args...) = flux
  return FluxBC{boundary}(flux_function)
end

function ValueBC(boundary, value::Number)
  @inline value_function(args...) = value
  return ValueBC{boundary}(value_function)
end

FluxBC(boundary, value) = FluxBC{boundary}(value)
ValueBC(boundary, value) = ValueBC{boundary}(value)
FluxBC(flux) = FluxBC(Unset, flux)
ValueBC(value) = ValueBC(Unset, value)

#
# Boundary condition API
#

"""
    FieldboundaryConditions(; top=TopBC, bottom=BottomBC)

Create an instance of `FieldBoundaryConditions` with `top` and `bottom` 
boundary conditions.
"""
function FieldBoundaryConditions(;
    bottom = FluxBC{Bottom}(zero_function),
       top = FluxBC{Top}(zero_function),
  )
  FieldBoundaryConditions(bottom, top)
end

"Returns `FieldBoundaryConditions` with zero flux at `top` and `bottom`."
ZeroFlux() = FieldBoundaryConditions() # the default

"Set the bottom boundary condition for `fld` in `model`."
function set_bottom_bc!(model, fld, bc) 
  field_bcs = getproperty(model.bcs, fld)
  field_bcs.bottom = bc
  return nothing
end

"Set the top boundary condition for `fld` in `model`."
function set_top_bc!(model, fld, bc) 
  field_bcs = getproperty(model.bcs, fld)
  field_bcs.top = bc
  return nothing
end

"""
    set_bc!(model, fld, bc)
    set_bc!(Boundary, model, fld, bc)

Add `bc` to `model` as (a) boundary condition(s) for `fld`.
If `Boundary = Top` or `Bottom` is specified, the boundary condition
is applied to the top or bottom accordingly.
"""
set_bc!(model, fld, bc::BoundaryCondition{B}) where B <: Bottom = set_bottom_bc!(model, fld, bc)
set_bc!(model, fld, bc::BoundaryCondition{B}) where B <: Top = set_top_bc!(model, fld, bc)

function set_bc!(model, fld, bc::BoundaryCondition{B}) where B <: Unset 
  throw("The boundary on which to apply the boundary condition must be specified!")
end

set_bc!(::Type{B}, fld, bc) where B = set_bc!(model, fld, convert(BoundaryCondition{B}, bc))

function set_bc!(model, fld, bcs::FieldBoundaryConditions)
  set_bottom_bc!(model, fld, bcs.bottom)
  set_top_bc!(model, fld, bcs.top)
  return nothing
end

function set_bc!(model, fld, bcs::Tuple{BoundaryCondition{B},BoundaryCondition{B}}) where B
  throw("Two boundary conditions given for one boundary!")
end

function set_bc!(model, fld, bcs::Tuple{BoundaryCondition{B},BoundaryCondition{B}}) where {B<:Unset}
  # Heuristically assume that first boundary condition is bottom and second is top
  set_bottom_bc!(model, fld, bcs[1])
  set_top_bc!(model, fld, bcs[2])
  return nothing
end

function set_bc!(model, fld, bcs::Tuple{BoundaryCondition{B},BoundaryCondition{C}}) where {B,C}
  # Assume boundary values have been determined
  set_bc!(model, fld, bcs[1])
  set_bc!(model, fld, bcs[2])
  return nothing
end



function set_flux_bc!(model, fld, bcs::Tuple)
  set_bottom_bc!(model, fld, FluxBC(Bottom, bcs[1]))
  set_top_bc!(model, fld, FluxBC(Top, bcs[2]))
  return nothing
end

function set_value_bc!(model, fld, bcs::Tuple)
  set_bottom_bc!(model, fld, ValueBC(Bottom, bcs[1]))
  set_top_bc!(model, fld, ValueBC(Top, bcs[2]))
  return nothing
end

"""
    set_bcs!(model; kwargs...)

Set boundary conditions for `model`.
"""
function set_bcs!(model; kwargs...)
  for (k, v) in kwargs
    set_bc!(model, k, v)
  end
  return nothing
end

function set_bcs!(model, fld, bcs...)
  for bc in bcs
    set_bc!(model, fld, bc)
  end
  return nothing
end

"""
    set_flux_bcs!(model; kwargs...)

Set flux boundary conditions for `model`.
"""
function set_flux_bcs!(model; kwargs...)
  for (k, v) in kwargs
    set_flux_bc!(model, k, v)
  end
  return nothing
end

"""
    set_value_bcs!(model; kwargs...)

Set value boundary conditions for `model`.
"""
function set_value_bcs!(model; kwargs...)
  for (k, v) in kwargs
    set_value_bc!(model, k, v)
  end
  return nothing
end
