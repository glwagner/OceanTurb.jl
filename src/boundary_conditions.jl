# Boundary conditions for OceanTurb.jl

import Base: convert, setproperty!

abstract type BCType end
struct Flux <: BCType end
struct Gradient <: BCType end
struct Value <: BCType end

"Return c₁, where c₁ = c₀ + ∂c/∂z * (z₁ - z₀) = c₀ + ∂c/∂z * Δz."
@inline linear_extrap(c₀, ∂c∂z, Δz) = c₀ + ∂c∂z * Δz

mutable struct BoundaryCondition{C, T}
    condition :: T
end

function Base.convert(::Type{BoundaryCondition{C, T1}}, bc::BoundaryCondition{C}) where {C, T1}
    return BoundaryCondition{C, T1}(convert(T1, bc.condition))
end

BoundaryCondition(bctype, bc) = BoundaryCondition{bctype, typeof(bc)}(bc)

GradientBoundaryCondition(bc) = BoundaryCondition(Gradient, bc)
FluxBoundaryCondition(bc) = BoundaryCondition(Flux, bc)
ValueBoundaryCondition(bc) = BoundaryCondition(Value, bc)

const BC = BoundaryCondition

@inline getbc(model, bc::BC{C, <:Number}) where C = bc.condition
@inline getbc(model, bc::BC{C, <:Function}) where C = bc.condition(model)

@inline gradient(bc::BC{<:Gradient}, model, args...) = getbc(model, bc)
@inline gradient(bc::BC{<:Flux},     model, κ, args...) = -getbc(model, bc) / κ
@inline gradient(bc::BC{<:Value},    model, κ, cᴺ, Δf, args...) = 2 * (getbc(model, bc) - cᴺ) / Δf

"""
    fill_bottom_ghost_cell!(c, κ, model, bc)

Update the bottom ghost cell of c given the boundary condition `bc`, `model`, and
diffusivity `kappa`. `kappa` is used only if a flux boundary condition
is specified.
"""
@inline function fill_bottom_ghost_cell!(bc, c, κ, model)
    @inbounds c[0] = linear_extrap(c[1], gradient(bc, model, κ, c[1], -Δf(c, 1)), -Δc(c, 1))
    return nothing
end

"""
    fill_top_ghost_cell!(c, κ, model, bc)

Update the top ghost cell of `c` given boundary condition `bc`, `model`, and
diffusivity `kappa`
"""
@inline function fill_top_ghost_cell!(bc, c, κ, model, N)
    @inbounds c[N+1] = linear_extrap(c[N], gradient(bc, model, κ, c[N], Δf(c, N)), Δc(c, N+1))
    return nothing
end

"""
    FieldBoundaryConditions(; top=TopBC, bottom=BottomBC)

Create an instance of `FieldBoundaryConditions` with `top` and `bottom`
"""
mutable struct FieldBoundaryConditions{B, T}
    bottom :: B
    top    :: T
end

"""
    ZeroFluxBoundaryConditions(T=Float64)

Construct `FieldBoundaryConditions` with a zero `FluxBoundaryCondition`
at top and bottom.
"""
function ZeroFluxBoundaryConditions(T=Float64)
    FieldBoundaryConditions(
        FluxBoundaryCondition(-zero(T)),
        FluxBoundaryCondition(-zero(T))
    )
end

"""
    DefaultBoundaryConditions(T=Float64)

Returns default oceanic boundary conditions: a zero
`GradientBoundaryCondition` on bottom and a zero
`FluxBoundaryCondition` on top.
"""
function DefaultBoundaryConditions(T=Float64)
    return FieldBoundaryConditions(
                GradientBoundaryCondition(-zero(T)),
                FluxBoundaryCondition(-zero(T))
                )
end

"""
    BoundaryConditions([T=Float64;] bottom = GradientBoundaryCondition(-zero(T)),
                                       top = FluxBoundaryCondition(-zero(T)))

Returns `FieldBoundaryConditions` with a `bottom` and `top` boundary condition.
The type `T` is only relevant for the default values of `bottom` and `top`.
"""
function BoundaryConditions(T=Float64; bottom = GradientBoundaryCondition(-zero(T)),
                                          top = FluxBoundaryCondition(-zero(T))
                           )
    return FieldBoundaryConditions(bottom, top)
end

"""
    set_bcs!(model; bcspecs...)

Set boundary conditions of model solution fields.
The keyword argument name must be the name of a model solution
and its value is a (bottom_bc, top_bc) tuple.

Example
========

julia> set_bcs!(model, c=(FluxBoundaryCondition(-1),
                          FluxBoundaryCondition(0))
                )
"""
function set_bcs!(model; bcspecs...)
    bcs = model.bcs
    for (ϕsym, bcs) in bcspecs
        ϕbcs = getproperty(model.bcs, ϕsym)
        setproperty!(ϕbcs, :bottom, bcs[1])
        setproperty!(ϕbcs, :top, bcs[2])
    end
    return nothing
end
