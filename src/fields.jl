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
import Base: +, *, -, setindex!, getindex, eachindex, lastindex, similar, eltype, length,
             @propagate_inbounds

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

Return a `Field{Face}` on `grid` with its data initialized to 0.
"""
function FaceField(A::DataType, grid)
    data = convert(A, fill(0, face_size(grid)))
    FaceField{typeof(data), typeof(grid), eltype(data)}(data, grid)
end

"""
    CellField(grid)

Return a `Field{Cell}` on `grid` with its data initialized to 0.
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

Return a `Field{Cell}` with its `data` located on the `grid`.
if `data` is an array, it must be broadcastable to `c.data`, where
`c` is a `Field{Cell}`.
"""
function CellField(data, grid)
    c = CellField(grid)
    set!(c, data)
    return c
end

"""
    FaceField(data, grid)

Return a `Field{Face}` with its `data` located on the `grid`.
if `data` is an array, it must be broadcastable to `f.data`, where
`f` is a `Field{Face}`.
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

function integral(c::CellField)
    total = 0
    for i in eachindex(c)
        @inbounds total += c[i] * Δf(c.grid, i)
    end
    return total
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

#
# Differential operators and such for fields
#

Δc(c::AbstractField, i_face) = Δc(c.grid, i_face)
Δf(c::AbstractField, i_cell) = Δf(c.grid, i_cell)

"""
    ∂z(a, i)

Return the discrete derivative of `a` at grid point `i`.

The derivative of a `Field{Cell}` is computed at face points,
and the derviative of a `Field{Face}` is computed at cell points.
"""
∂z(a, i) = throw("∂z is not defined for arbitrary fields.")

"Return ∂c/∂z at face index i."
@propagate_inbounds ∂z(c::CellField, i) = (c.data[i] - c.data[i-1]) / Δc(c.grid, i)

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
# A bunch of (unsafe) diffusive flux operators
#

# ∇K∇c for c::CellField
K∂z(K, c, i) = @inbounds K*∂z(c, i)
∇K∇c(Kᵢ₊₁, Kᵢ, c, i)            = @inbounds ( K∂z(Kᵢ₊₁, c, i+1) -    K∂z(Kᵢ, c, i)     ) /    Δf(c, i)
∇K∇c_top(Kᴺ, c, top_flux)       = @inbounds (     -top_flux     - K∂z(Kᴺ, c, c.grid.N) ) / Δf(c, c.grid.N)
∇K∇c_bottom(K₂, c, bottom_flux) = @inbounds (   K∂z(K₂, c, 2)   +     bottom_flux      ) /    Δf(c, 1)

## Top and bottom flux estimates for constant (Dirichlet) boundary conditions
bottom_flux(K, c, c_bndry, Δf) = -2K*( bottom(c) - c_bndry ) / Δf # -K*∂c/∂z at the bottom
top_flux(K, c, c_bndry, Δf)    = -2K*(  c_bndry  -  top(c) ) / Δf # -K*∂c/∂z at the top

∇K∇c_top(Kᴺ⁺¹, Kᴺ, c, bc, model) = ∇K∇c_top(Kᴺ, c, -Kᴺ⁺¹*getbc(model, bc))
∇K∇c_bottom(K₂, K₁, c, bc, model) = ∇K∇c_bottom(K₂, c, -K₁*getbc(model, bc))

"Return the total flux (advective + diffusive) across face i."
@propagate_inbounds flux(w, κ, c, i) = w * onface(c, i) - κ * ∂z(c, i)
top_flux_div(wtop, κtop, c) = @inbounds -flux(wtop, κtop, c, c.grid.N) / Δf(c, c.grid.N)
bottom_flux_div(wbottom, κbottom, c) = @inbounds flux(wbottom, κbottom, c, 1) / Δf(c, 1)

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
@propagate_inbounds onface(c::CellField, i) = 0.5*(c.data[i] + c.data[i-1])
@propagate_inbounds onface(f::FaceField, i) = f[i]

"""
    oncell(f, i)

Return the interpolation of `f` onto cell point `i`.
"""
@propagate_inbounds oncell(f::FaceField, i) = 0.5*(f.data[i+1] + f.data[i])
@propagate_inbounds oncell(c::CellField, i) = c[i]
