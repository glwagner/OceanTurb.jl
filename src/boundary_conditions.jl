# Boundary conditions for OceanTurb.jl

import Base: convert

abstract type BCType end
struct Flux <: BCType end
struct Gradient <: BCType end
struct Value <: BCType end

#abstract type BoundaryCondition{B<:Boundary} end

@inline zero_function(args...) = 0

struct BoundaryCondition{C<:BCType, T}
    condition :: T
end

const BC = BoundaryCondition

BoundaryCondition(C, bc) = BoundaryCondition{C, typeof(bc)}(bc)

getbc(model, bc::BC{C, <:Function}) where C = bc.condition(model)
getbc(model, bc::BC{C, <:Number}) where C = bc.condition

"""
    FluxBoundaryCondition(boundary, flux)

Constuct a boundary condition that specifies the flux
of some field on a boundary. If `flux` is a function,
its arguments must be synced with the expection of `Model`.
"""
FluxBoundaryCondition(bc) = BoundaryCondition(Flux, bc)

"""
    ValueBoundaryCondition(boundary, flux)

Constuct a boundary condition that specifies the value
of some field on a boundary. If `flux` is a function,
its arguments must be synced with the expection of `Model`.
"""
ValueBoundaryCondition(bc) = BoundaryCondition(Value, bc)

"""
    GradientBoundaryCondition(boundary, flux)

Constuct a boundary condition that specifies the gradient
of some field on a boundary. If `flux` is a function,
its arguments must be synced with the expection of `Model`.
"""
GradientBoundaryCondition(bc) = BoundaryCondition(Gradient, bc)


mutable struct FieldBoundaryConditions
    bottom::BoundaryCondition
    top::BoundaryCondition
end

#
# Boundary condition API
#

"""
    FieldBoundaryConditions(; top=TopBC, bottom=BottomBC)

Create an instance of `FieldBoundaryConditions` with `top` and `bottom`
boundary conditions.
"""
function FieldBoundaryConditions(;
    bottom = FluxBoundaryCondition(0),
    top = FluxBoundaryCondition(0),
    )
    FieldBoundaryConditions(bottom, top)
end

"Returns `FieldBoundaryConditions` with zero flux at `top` and `bottom`."
ZeroFluxBoundaryConditions() = FieldBoundaryConditions() # the default

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

function set_bc!(model, fld, bcs::FieldBoundaryConditions)
  set_bottom_bc!(model, fld, bcs.bottom)
  set_top_bc!(model, fld, bcs.top)
  return nothing
end

function set_top_flux_bc!(model, fld, bc)
    set_top_bc!(model, fld, FluxBoundaryCondition(bc))
    return nothing
end

function set_bottom_flux_bc!(model, fld, bc)
    set_bottom_bc!(model, fld, FluxBoundaryCondition(bc))
    return nothing
end

"""
    set_bcs!(model, fld, bcs::Tuple)

Set boundary conditions on `fld` as  `bcs = (bottom_bc, top_bc)`

    set_bcs!(model; kwargs...)

Set boundary conditions, where the keywords are fields of `model.solution` and
their arguments are tuples of (bottom, top) boundary condition.
"""
function set_bcs!(model, fld, bcs::Tuple)
  # First boundary condition is bottom and second is top
  set_bottom_bc!(model, fld, bcs[1])
  set_top_bc!(model, fld, bcs[2])
  return nothing
end

function set_bcs!(model; kwargs...)
    for (k, v) in kwargs
        set_bcs!(model, k, v)
    end
    return nothing
end

"""
    set_flux_bcs!(model, fld, bcs::Tuple)

Set flux boundary conditions on `fld` as  `bcs = (bottom_bc, top_bc)`

    set_flux_bcs!(model; kwargs...)

Set flux boundary conditions for `model`.
"""
function set_flux_bcs!(model; kwargs...)
    for (k, v) in kwargs
        set_flux_bcs!(model, k, v)
    end
    return nothing
end

function set_flux_bcs!(model, fld, bcs)
    set_top_flux_bc!(model, fld, bc[1])
    set_bottom_flux_bc!(model, fld, bc[2])
    return nothing
end
