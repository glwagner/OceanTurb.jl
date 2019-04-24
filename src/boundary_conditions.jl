# Boundary conditions for OceanTurb.jl

import Base: convert, setproperty!

abstract type BCType end
struct Flux <: BCType end
struct Gradient <: BCType end
struct Value <: BCType end

#abstract type BoundaryCondition{B<:Boundary} end

@inline zero_function(args...) = 0

"Return c₁, where c₁ = c₀ + ∂c/∂z * (z₁ - z₀) = c₀ + ∂c/∂z * Δz."
linear_extrapolation(c₀, ∂c∂z, Δz) = c₀ + ∂c∂z * Δz

getbc(model, bc::Function) = bc(model)
getbc(model, bc::Number) = bc

struct BoundaryCondition
    gradient :: Function
end

function FluxBoundaryCondition(bc)
    gradient(c, Δf, κ, model) = -getbc(model, bc) / κ
    return BoundaryCondition(gradient)
end

function GradientBoundaryCondition(bc)
    gradient(c, Δf, κ, model) = getbc(model, bc)
    return BoundaryCondition(gradient)
end

function ValueBoundaryCondition(bc)
    gradient(cᴺ, Δf, κ, model) = 2 * (getbc(model, bc) - cᴺ) / Δf
    return BoundaryCondition(gradient)
end

"""
    fill_bottom_ghost_cell!(c, κ, model, bc)

Update the bottom ghost cell of c given the boundary condition `bc`, `model`, and
diffusivity `kappa`. `kappa` is used only if a flux boundary condition
is specified.
"""
function fill_bottom_ghost_cell!(c, κ, model, bc::BoundaryCondition)
    @inbounds begin
        c[0] = linear_extrapolation(
            c[1], bc.gradient(c[1], -Δf(c.grid, 1), κ, model), -Δc(c.grid, 1))
    end
    return nothing
end

"""
    fill_top_ghost_cell!(c, κ, model, bc)

Update the top ghost cell of `c` given boundary condition `bc`, `model`, and
diffusivity `kappa`
"""
function fill_top_ghost_cell!(c, κ, model, bc::BoundaryCondition)
    N = c.grid.N
    @inbounds begin
        c[N+1] = linear_extrapolation(
            c[N], bc.gradient(c[N], Δf(c.grid, N), κ, model), Δc(c.grid, N+1))
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

mutable struct FieldBoundaryConditions
    bottom :: BoundaryCondition
    top    :: BoundaryCondition
end

#
# Boundary condition API
#

"""
    FieldBoundaryConditions(; top=TopBC, bottom=BottomBC)

Create an instance of `FieldBoundaryConditions` with `top` and `bottom`
boundary conditions.
"""
function FieldBoundaryConditions(; bottom = FluxBoundaryCondition(0),
                                      top = FluxBoundaryCondition(0))
    FieldBoundaryConditions(bottom, top)
end

"Returns `FieldBoundaryConditions` with zero flux at `top` and `bottom`."
ZeroFluxBoundaryConditions() = FieldBoundaryConditions() # the default
