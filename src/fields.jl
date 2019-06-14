#= Fields for OceanTurb.jl

A `Field` is an abstraction of a variable or function defined on a staggered grid.
Our hope is that it simplifies differential operations on a staggered grid.

OceanTurb.jl solves one-dimensional PDEs on a staggered grid.
The geometry of a grid with `N=3` is

```
      ▲ z
      |

         i=4           *
                j=4   ===       ▲
         i=3           *        | Δf (i=3)
                j=3   ---       ▼
         i=2           *    ▲
                j=2   ---   | Δc (j=2)
         i=1           *    ▼
                j=1   ===
         i=0           *
```

where the i's index cells and the j's index faces.
The variable Δc gives the separation between
cell centers, and Δf gives the separation between faces.
Ghost cells at i=0 and i=N+1 bound the domain.

There are two types of fields:

  1. Fields defined at cell centers with dimension `N+2`: `Field{Cell}`
  2. Fields defined at cell interfaces with dimension `N+1`: `Field{Face}`
=#
import Base: +, *, -, ^, setindex!, getindex, eachindex, lastindex, similar,
             eltype, length, @propagate_inbounds

import Statistics: mean

default_arraytype(T) = Array{T, 1}

struct CellField{A, G, T} <: AbstractField{A, G, T}
    data :: OffsetArray{T, 1, A}
    grid :: G
end

struct FaceField{A, G, T} <: AbstractField{A, G, T}
    data :: A
    grid :: G
end

arraytype(::AbstractField{A}) where A = A
eltype(::AbstractField{A}) where A = eltype(A)

#
# (legacy) Field Location and Constructors
#

abstract type FieldLocation end
struct Cell <: FieldLocation end
struct Face <: FieldLocation end

Field(::Type{Face}, grid) = FaceField(grid)
Field(::Type{Cell}, grid) = CellField(grid)

"""
    FaceField(grid)

Return a `FaceField` on `grid` with its data initialized to 0.
"""
function FaceField(A::DataType, grid)
    data = convert(A, fill(0, face_size(grid)))
    FaceField{typeof(data), typeof(grid), eltype(data)}(data, grid)
end

"""
    CellField(grid)

Return a `CellField` on `grid` with its data initialized to 0.
"""
function CellField(A::DataType, grid)
    data = convert(A, fill(0, cell_size(grid)))
    offset_data = OffsetArray(data, 0:grid.N+1)
    CellField{typeof(data), typeof(grid), eltype(data)}(offset_data, grid)
end

CellField(grid) = CellField(default_arraytype(eltype(grid)), grid)
FaceField(grid) = FaceField(default_arraytype(eltype(grid)), grid)

"""
    CellField(data, grid)

Return a `CellField` with its `data` located on the `grid`.
if `data` is an array, it must be broadcastable to `c.data`, where
`c` is a `CellField`.
"""
function CellField(data, grid)
    c = CellField(grid)
    set!(c, data)
    return c
end

"""
    FaceField(data, grid)

Return a `FaceField` with its `data` located on the `grid`.
if `data` is an array, it must be broadcastable to `f.data`, where
`f` is a `FaceField`.
"""
function FaceField(data, grid)
    f = FaceField(grid)
    set!(f, data)
    return f
end

#
# Basic 'Field' functionality
#

data(c::FaceField) = c.data
data(c::CellField) = view(c.data, 1:c.grid.N)

nodes(c::CellField) = c.grid.zc
nodes(f::FaceField) = f.grid.zf

length(c::CellField) = c.grid.N
length(f::FaceField) = f.grid.N + 1

height(c::AbstractField) = height(c.grid)

# All indices
eachindex(c::CellField) = 1:c.grid.N
eachindex(f::FaceField) = 1:f.grid.N + 1

lastindex(c::CellField) = c.grid.N
lastindex(f::FaceField) = f.grid.N + 1

# Interior indices, omitting boundary-adjacent values
interiorindices(c::CellField) = 2:c.grid.N - 1
interiorindices(f::FaceField) = 2:f.grid.N
boundaryindices(c::CellField) = (1, c.grid.N)

# Sugary sweet: access indices of c.data by indexing into c.
@propagate_inbounds getindex(c::AbstractField, inds...) = getindex(c.data, inds...)
@propagate_inbounds setindex!(c::AbstractField, d, inds...) = setindex!(c.data, d, inds...)
@propagate_inbounds setindex!(c::AbstractField, d::AbstractField, inds...) = setindex!(c.data, d.data, inds...)

#
# Ways to specify a field's data
#

set!(c::AbstractField, data::Number) = fill!(c.data, data)
set!(c::AbstractField{Ac, G}, d::AbstractField{Ad, G}) where {Ac, Ad, G} = c.data .= convert(Ac, d.data)
set!(c::FaceField, fcn::Function) = c.data .= fcn.(nodes(c))

function set!(c::CellField, func::Function)
    data = func.(nodes(c))
    set!(c, data)
    # Set ghost points to get approximation to first derivative at boundary
    data_bottom = func(c.grid.zf[1])
    data_top = func(c.grid.zf[end])

    # Set ghost values so that
    # ∂z(c, 1) = (c[1] - c[0]) / Δc(c, 1) = (c[1] - c_bottom) / 0.5*Δc(c, 1)
    #
    # and
    # ∂z(c, N+1) = (c[N+1] - c[N]) / Δc(c, N+1) = (c_top - c[N]) / 0.5*Δc(c, N)

    N = c.grid.N
    @inbounds begin
        c[0] = c[1] - 2 * (c[1] - data_bottom)
        c[N+1] = c[N] + 2 * (data_top - c[N])
    end

    return nothing
end

set!(f::FaceField, data::AbstractArray) = f.data .= data

function set!(c::CellField, data::AbstractArray)
    for i in eachindex(data)
        @inbounds c[i] = data[i]
    end
    # Default boundary conditions if data is not an OffsetArray
    typeof(data) <: OffsetArray || set_default_bcs!(c)
    return nothing
end

function set_default_bcs!(c)
    @inbounds begin
        c[0] = c[1]
        c[c.grid.N+1] = c[c.grid.N]
    end
    return nothing
end

function integral(fn::Function, c::CellField)
    total = zero(eltype(c))
    for i in eachindex(c)
        @inbounds total += fn(c[i]) * Δf(c, i)
    end
    return total
end

integral(c::CellField) = integral(x->x, c)

"""
    mean([f], c::CellField)

Compute the mean of the field `c` over its domain,
applying the function `f` to each element.
`f` is the identity function f(x) = x by default.
"""
mean(fn::Function, c::CellField) = integral(fn, c) / height(c)
mean(c::CellField) = mean(x->x, c)


function integrate_range(c::CellField, i₁::Int, i₂::Int)
    total = 0
    for i = i₁:i₂
        @inbounds total += c[i] * Δf(c.grid, i)
    end
    return total
end

function integral(c::CellField, z₋, z₊=0)

    @assert z₊ > c.grid.zf[1] "Integration region lies outside the domain."
    @assert z₊ > z₋ "Invalid integration range: upper limit greater than lower limit."

    # Find region bounded by the face ≤ z₊ and the face ≤ z₁
    i₁ = searchsortedfirst(c.grid.zf, z₋) - 1
    i₂ = searchsortedfirst(c.grid.zf, z₊) - 1

    if i₂ ≠ i₁
        # Calculate interior integral, recalling that the
        # top interior cell has index i₂-2.
        total = integrate_range(c, i₁+1, i₂-1)

        # Add contribution to integral from fractional bottom part,
        # if that region is a part of the grid.
        if i₁ > 0
            total += c[i₁] * (c.grid.zf[i₁+1] - z₋)
        end

        # Add contribution to integral from fractional top part
        total += c[i₂] * (z₊ - c.grid.zf[i₂])
    else
        total = c[i₁] * (z₊ - z₋)
    end

    return total
end

function set!(c1::CellField{A1, G2}, c2::CellField{A2, G2}) where {A1, G1, A2, G2}
    @assert height(c1) == height(c2) "Physical domains differ between the two fields."
    for i in eachindex(c1)
        @inbounds c1[i] = integral(c2, c1.grid.zf[i], c1.grid.zf[i+1]) / Δf(c1, i)
    end
    return nothing
end

similar(c::CellField) = CellField(c.grid)
similar(f::FaceField) = FaceField(f.grid)

# Define +, -, and * on fields as element-wise calculations on their data. This
# is only true for fields of the same type. So far, we haven't found use for
# these sweets because we tend to write element-wise kernels for operations.
for op in (:+, :-, :*)
    @eval begin
        # +, -, * a Field by a Number on the left
        function $op(num::Number, f::AbstractField)
            ff = similar(f)
            @. ff.data = $op(num, f.data)
            ff
        end

        # +, -, * a Field by a Number on the right.
        $op(f::AbstractField, num::Number) = $op(num, f)

        # Binary two-field operations
        function $op(f1::F, f2::F) where {F <: AbstractField}
            f3 = similar(f1)
            @. f3.data = $op(f1.data, f2.data)
            f3
        end
    end
end

function ^(c::AbstractField, b::Number)
    d = similar(c)
    set!(d, c.data.^b)
    return d
end

#
# Differential operators and such for fields
#

@propagate_inbounds Δc(c::AbstractField, i_face) = Δc(c.grid, i_face)
@propagate_inbounds Δf(c::AbstractField, i_cell) = Δf(c.grid, i_cell)

"""
    ∂z(a, i)

Return the discrete derivative of `a` at grid point `i`.

The derivative of a `CellField` is computed at face points,
and the derviative of a `FaceField` is computed at cell points.
"""
∂z(a, i) = throw("∂z is not defined for arbitrary fields.")

"Return ∂c/∂z at face index i."
@propagate_inbounds ∂z(c::CellField, i) = (c.data[i] - c.data[i-1]) / Δc(c, i)

"Return ∂c/∂z at face index i."
@propagate_inbounds ∂z(c::FaceField, i) = (c.data[i+1] - c.data[i]) / Δc(c, i)
@propagate_inbounds ∂²z(c::AbstractField, i) = (∂z(c, i+1) - ∂z(c, i)) / Δf(c, i)

"Calculate `f = ∂c/∂z` in the grid interior."
function ∂z!(f::FaceField, c::CellField)
    for i = eachindex(f)
        @inbounds f.data[i] = ∂z(c, i)
    end
    return nothing
end

"Calculate `c = ∂f/∂z` in the grid interior."
function ∂z!(c::CellField, f::FaceField)
    for i = eachindex(c)
        @inbounds c.data[i] = ∂z(f, i)
    end
    return nothing
end

"Return the `FaceField` ∂c/∂z, where `c` is a `CellField`."
function ∂z(c::CellField)
    f = FaceField(c.grid)
    ∂z!(f, c)
    return f
end

"Return the `CellField` ∂f/∂z, where `f` is a `FaceField`."
∂z(a) = throw("∂z is not defined for arbitrary fields.")

function ∂z(f::FaceField)
    c = CellField(f.grid)
    ∂z!(c, f)
    return c
end

#
# Advection and diffusion operators
#

@propagate_inbounds K∂z(K, ϕ, i) = K * ∂z(ϕ, i)

"Return the diffusive flux divergence at cell i."
@propagate_inbounds ∇K∇ϕ(Kᵢ₊₁, Kᵢ, ϕ, i) = (K∂z(Kᵢ₊₁, ϕ, i+1) - K∂z(Kᵢ, ϕ, i)) / Δf(ϕ, i)

"Return the upwind advective flux divergence at cell i for M<0."
@propagate_inbounds ∂zM(Mᵢ₊₁, Mᵢ, ϕ, i) = (Mᵢ₊₁ * ϕ[i+1] - Mᵢ * ϕ[i]) / Δc(ϕ, i+1)

"Return the total flux (advective + diffusive) across face i."
@propagate_inbounds flux(M, K, ϕ, i) = M * onface(ϕ, i) - K * ∂z(ϕ, i)

#
# Convenience functions
#

top(a) = a
top(a::AbstractArray) = @inbounds a[end]
top(a::CellField) = @inbounds a[a.grid.N]
top(a::FaceField) = @inbounds a[a.grid.N+1]

bottom(a) = throw("bottom(a) Not implemented for typeof(a) = $(typeof(a)).")
bottom(a::Number) = a
bottom(a::Union{AbstractField, AbstractArray}) = @inbounds a[1]

"""
    onface(c, i)

Return the interpolation of `c` onto face point `i`.
"""
@propagate_inbounds onface(c::CellField, i) = (c.data[i] + c.data[i-1])/2
@propagate_inbounds onface(f::FaceField, i) = f[i]

"""
    oncell(f, i)

Return the interpolation of `f` onto cell point `i`.
"""
@propagate_inbounds oncell(f::FaceField, i) = (f.data[i+1] + f.data[i])/2
@propagate_inbounds oncell(c::CellField, i) = c[i]

"""
    absolute_error(c, d, p=2)

Compute the absolute error between `c` and `d` with norm `p`, defined as

``\\mathrm{abs \\, error} = \\left ( L^{-1} \\int_{-L}^0 (c-d)^p \\, \\mathrm{d} z \\right )^(1/p)``.
"""
function absolute_error(c::CellField, d::CellField, p=2)
    if length(c) != length(d)
        𝒹 = similar(c)
        set!(𝒹, d)
    else
        𝒹 = d
    end

    total = zero(eltype(c))
    for i in eachindex(c)
        @inbounds total += (c[i] - 𝒹[i])^p * Δf(c, i)
    end

    return  ( total / height(c) )^(1/p)
end

"""
    relative_error(c, d, p=2)

Compute the relative error between `c` and `d` with norm `p`, defined as

```math
\\beq
\\mathrm{rel \\, error} = \\frac{\\left ( int_{-L}^0 (c-d)^p \\, \\mathrm{d} z \\right )^(1/p)}
                             {\\left ( int_{-L}^0 d^p \\, \\mathrm{d} z \\right )^(1/p)}
\\eeq
```
"""
relative_error(c::CellField, d::CellField, p=2) = absolute_error(c, d, p) / mean(x -> x^p, d)^(1/p)
