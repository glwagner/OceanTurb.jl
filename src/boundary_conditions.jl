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

Construct a boundary condition that specifies a field's flux
of some field on a boundary.
"""
FluxBoundaryCondition(bc) = BoundaryCondition(Flux, bc)

"""
    ValueBoundaryCondition(boundary, value)

Construct a boundary condition that specifies a field's value
of some field on a boundary.
"""
ValueBoundaryCondition(bc) = BoundaryCondition(Value, bc)

"""
    GradientBoundaryCondition(boundary, gradient)

Construct a boundary condition that specifies a field's gradient
on a boundary.
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

"Return c₁, where c₁ = c₀ + ∂c/∂z * (z₁ - z₀) = c₀ + ∂c/∂z * Δz."
linear_extrapolation(c₀, ∂c∂z, Δz) = c₀ + ∂c∂z * Δz

get_gradient(c₀, Δf, κ, model, bc::BoundaryCondition{<:Gradient}) = getbc(model, bc)
get_gradient(c₀, Δf, κ, model, bc::BoundaryCondition{<:Flux}) = -getbc(model, bc) / κ

function get_gradient(cᴺ, Δf, κ, model, bc::BoundaryCondition{<:Value})
    return 2 * (getbc(model, bc) - cᴺ) / Δf
end

"""
    fill_top_ghost_cell!(c, κ, model, bc)

Update the top ghost cell of `c` given boundary condition `bc`, `model`, and
diffusivity `kappa`
"""
function fill_top_ghost_cell!(c, κ, model, bc)
    @inbounds begin
        c[c.grid.N+1] = linear_extrapolation(
            c[c.grid.N],
            get_gradient(c[c.grid.N], Δf(c.grid, c.grid.N), κ, model, bc),
            Δc(c.grid, c.grid.N+1))
    end
    return nothing
end

"""
    fill_bottom_ghost_cell!(c, κ, model, bc)

Update the bottom ghost cell of c given the boundary condition `bc`, `model`, and
diffusivity `kappa`. `kappa` is used only if a flux boundary condition
is specified.
"""
function fill_bottom_ghost_cell!(c, κ, model, bc)
    @inbounds begin
        c[0] = linear_extrapolation(c[1],
                                    get_gradient(c[1], -Δf(c.grid, 1), κ, model, bc),
                                    -Δc(c.grid, 1))
    end
    return nothing
end

"""
    fill_ghost_cells!(c, κ, model, bc)

Update the top and bottom ghost cells of `c` for `model`.
"""
function fill_ghost_cells!(c, κbottom, κtop, model, fieldbcs)
    fill_bottom_ghost_cell!(c, κbottom, model, fieldbcs.bottom)
    fill_top_ghost_cell!(c, κtop, model, fieldbcs.top)
    return nothing
end

# Convenience methods for pre-computed fluxes
get_gradient(κ, flux) = - flux / κ

function fill_top_ghost_cell_flux!(c, κ, flux)
    @inbounds begin
        c[c.grid.N+1] = linear_extrapolation(c[c.grid.N], get_gradient(κ, flux),
                                             Δc(c.grid, c.grid.N+1))
    end
    return nothing
end
