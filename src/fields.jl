#= Fields for OceanTurb.jl

A `Field` is an abstraction of a variable or function defined on a staggered grid.
Our hope is that it simplifies differential operations on a staggered grid.

OceanTurb.jl solves one-dimensional PDEs on a staggered grid.
The geometry of a grid with `nz=3` is 

```
      ▲ z 
      |   
        
                j=4   ===       ▲              
         i=3           *        | dzf (i=3)
                j=3   ---       ▼
         i=2           *    ▲            
                j=2   ---   | dzc (j=2) 
         i=1           *    ▼  
                j=1   ===     
```

where the i's index cells and the j's index faces. 
The variable dzc gives the separation between
cell centers, and dzf gives the separation between faces.

There are two types of fields:

  1. Fields defined at cell centers with dimension `nz`: `Field{Cell}`
  2. Fields defined at cell faces with dimension `nz+1`: `Field{Face}`

From the standpoint of designing new turbulence closures, 
the most important function that we output is `∂z`.
=#

import Base: +, *, -, setindex!, getindex, eachindex, similar, eltype, length

abstract type FieldLocation end
struct Cell <: FieldLocation end
struct Face <: FieldLocation end

struct Field{L,A,G} <: AbstractField{A,G}
  data::A
  grid::G
  function Field(Location, data, grid)
    new{Location,typeof(data),typeof(grid)}(data, grid)
  end
end

const CellField = Field{Cell}
const FaceField = Field{Face}

#
# Field Constructors
#

Field(::Type{Face}, grid) = FaceField(grid)
Field(::Type{Cell}, grid) = CellField(grid)

"""
    FaceField(grid)

Return a `Field{Face}` on `grid` with its data initialized to 0.
"""
function FaceField(grid)
  A = arraytype(grid)
  data = convert(A, fill(0, face_size(grid)))
  Field(Face, data, grid)
end

"""
    CellField(grid)

Return a `Field{Cell}` on `grid` with its data initialized to 0.
"""
function CellField(grid)
  A = arraytype(grid)
  data = convert(A, fill(0, cell_size(grid)))
  Field(Cell, data, grid)
end

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

zdata(c::CellField) = c.grid.zc
zdata(f::FaceField) = f.grid.zf

length(c::CellField) = cell_length(c.grid)
length(f::FaceField) = face_length(f.grid)

# All indices
eachindex(f::AbstractField) = eachindex(f.data)

# Interior indices, omitting boundary-adjacent values
interior(c::CellField) = 2:c.grid.nz-1 
interior(f::FaceField) = 2:f.grid.nz 

# Sugary sweet: access indices of c.data by indexing into c.
getindex(c::AbstractField, inds...) = getindex(c.data, inds...)
setindex!(c::AbstractField, d, inds...) = setindex!(c.data, d, inds...)
setindex!(c::AbstractField, d::Field, inds...) = setindex!(c.data, d.data, inds...)

set!(c::AbstractField, data::Number) = fill!(c.data, data)
set!(c::AbstractField{A}, data::AbstractArray) where A = c.data .= convert(A, data)
set!(c::AbstractField{A}, data::Function) where A = c.data .= convert(A, data.(zdata(c)))
set!(c::AbstractField{Ac,G}, d::AbstractField{Ad,G}) where {Ac,Ad,G} = c.data .= convert(Ac, d.data)

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
      function $op(f1::Field{L}, f2::Field{L}) where L
          f3 = similar(f1)
          @. f3.data = $op(f1.data, f2.data)
          f3  
      end 
  end 
end 

#
# Differential operators and such for fields
# 

"Return the cell spacing at index i."
dzc(c, i_face) = c.grid.dzc[i_face]
dzc(c::AbstractField{A,G}, i) where {A,G<:UniformGrid} = c.grid.dzc

"Return the face spacing at index i."
dzf(c, i_cell) = c.grid.dzf[i_cell]
dzf(c::AbstractField{A,G}, i) where {A,G<:UniformGrid} = c.grid.dzf

"""
    ∂z(a, i)

Return the discrete derivative of `a` at grid point `i`.

The derivative of a `Field{Cell}` is computed at face points,
and the derviative of a `Field{Face}` is computed at cell points.
"""
∂z(a, i) = throw("∂z is not defined for arbitrary fields.")

"Return ∂c/∂z at face index i."
∂z(c::CellField, i) = (c.data[i] - c.data[i-1]) / dzc(c, i)

"Return ∂c/∂z at face index i."
∂z(c::FaceField, i) = (c.data[i+1] - c.data[i]) / dzc(c, i)
∂²z(c::AbstractField, i) = (∂z(c, i+1) - ∂z(c, i)) / dzf(c, i)

"Calculate `f = ∂c/∂z` in the grid interior."
function ∂z!(f::FaceField, c::CellField)
  for i = interior(f)
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

# Convenience functions
top(a) = throw("top(a) Not implemented for typeof(a) = $(typeof(a)).")
top(a::Number) = a
top(a::AbstractArray) = a[end]
top(a::Field) = a.data[end]

bottom(a) = throw("bottom(a) Not implemented for typeof(a) = $(typeof(a)).")
bottom(a::Number) = a
bottom(a::AbstractArray) = a[1]
bottom(a::Field) = a.data[1]

"""
    avz(c, i)

Return the average of a at index i.
The average of a `Field{Cell}` is computed at face points.
The average of a `Field{Face}` is computed at cell points.
"""
avz(a, i) = throw("avz not defined for arbitrary fields.")
avz(f::FaceField, i) = 0.5*(f.data[i+1] + f.data[i])

function avz(c::CellField, i) 
  if i > 1
    return 0.5*(c.data[i] + c.data[i-1])
  else
    return c.data[1]
  end
end
