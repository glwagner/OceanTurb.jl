# Boundary conditions for OceanTurb.jl

import Base: convert, setproperty!

abstract type BCType end
struct Flux <: BCType end
struct Gradient <: BCType end
struct Value <: BCType end

"Return c₁, where c₁ = c₀ + ∂c/∂z * (z₁ - z₀) = c₀ + ∂c/∂z * Δz."
linear_extrap(c₀, ∂c∂z, Δz) = c₀ + ∂c∂z * Δz

abstract type AbstractBoundaryCondition end

struct BoundaryCondition{C, T} <: AbstractBoundaryCondition
    condition :: T
end

BoundaryCondition(bctype, bc) = BoundaryCondition{bctype, typeof(bc)}(bc)
GradientBoundaryCondition(bc) = BoundaryCondition(Gradient, bc)
FluxBoundaryCondition(bc) = BoundaryCondition(Flux, bc)
ValueBoundaryCondition(bc) = BoundaryCondition(Value, bc)

struct ZeroFluxBoundaryCondition <: AbstractBoundaryCondition end

const BC = BoundaryCondition
const ZFBC = ZeroFluxBoundaryCondition

getbc(model, ::ZFBC) = 0
getbc(model, bc::BC{C, <:Number}) where C = bc.condition
getbc(model, bc::BC{C, <:Function}) where C = bc.condition(model)

gradient(bc::BC{<:Gradient}, model, args...) = getbc(model, bc)
gradient(bc::BC{<:Flux},     model, κ, args...) = -getbc(model, bc) / κ
gradient(bc::BC{<:Value},    model, κ, cᴺ, Δf, args...) = 2 * (getbc(model, bc) - cᴺ) / Δf

"""
    fill_bottom_ghost_cell!(c, κ, model, bc)

Update the bottom ghost cell of c given the boundary condition `bc`, `model`, and
diffusivity `kappa`. `kappa` is used only if a flux boundary condition
is specified.
"""
function fill_bottom_ghost_cell!(bc, c, κ, model)
    @inbounds c[0] = linear_extrap(c[1], gradient(bc, model, κ, c[1], -Δf(c, 1)), -Δc(c, 1))
    return nothing
end

"""
    fill_top_ghost_cell!(c, κ, model, bc)

Update the top ghost cell of `c` given boundary condition `bc`, `model`, and
diffusivity `kappa`
"""
function fill_top_ghost_cell!(bc, c, κ, model, N)
    @inbounds c[N+1] = linear_extrap(c[N], gradient(bc, model, κ, c[N], Δf(c, N)), Δc(c, N+1))
    return nothing
end

fill_bottom_ghost_cell!(::ZFBC, c, args...) = @inbounds c[0] = c[1]
fill_top_ghost_cell!(::ZFBC, c, κ, model, N) = @inbounds c[N+1] = c[N]

mutable struct FieldBoundaryConditions
    bottom :: AbstractBoundaryCondition
    top    :: AbstractBoundaryCondition
end

"""
    FieldBoundaryConditions(; top=TopBC, bottom=BottomBC)

Create an instance of `FieldBoundaryConditions` with `top` and `bottom`
boundary conditions.
"""
function FieldBoundaryConditions(; bottom = ZeroFluxBoundaryCondition(),
                                      top = ZeroFluxBoundaryCondition())
    FieldBoundaryConditions(bottom, top)
end

ZeroFluxBoundaryConditions() = FieldBoundaryConditions()
