#= Fields for OceanTurb.jl

A `Field` is an abstraction of a variable or function defined on a staggered grid.
Our hope is that it simplifies differential operations on a staggered grid.

OceanTurb.jl solves one-dimensional PDEs on a staggered grid.
The geometry of a grid with nz=3 is 

```
      ^ z 
      |   
        
                j=4   ===       ^              
         i=3           *        | dzf (i=3)
                j=3   ---       v
         i=2           *    ^            
                j=2   ---   | dzc (j=2) 
         i=1           *    v  
                j=1   ===     
```
 

where the i's index cells and the j's index faces. 
The variable dzc gives the separation between
cell centers, and dzf gives the separation between faces.

There are two types of fields:

  1. `CellFields` defined at cell centers with dimension `nz`, 
  2. `FaceFields` defined at cell faces with dimension `nz+1`.

From the standpoint of designing new turbulence closures, 
the most important function that we output is `∂z`.
=#

import Base: +, *, -, setindex!, getindex, eachindex, similar, setproperty!, eltype

"""
    CellField(data, grid)

Return a `CellField` with `data` at the cell points `grid.zc`.
"""
struct CellField{A,G} <: Field{A,G} 
  data::A
  grid::G
end

"""
    FaceField(data, grid)

Return a `FaceField` with `data` at the face points `grid.zf`.
"""
struct FaceField{A,G} <: Field{A,G}
  data::A
  grid::G
end

#
# Constructors for Cell and Face fields
#

function CellField(data::Array, grid)
  A = arraytype(grid)
  CellField{A,typeof(grid)}(data, grid)
end

function FaceField(data::Array, grid::Grid{T,AG}) where {T,AG} 
  A = arraytype(grid)
  FaceField{A,typeof(grid)}(data, grid)
end

function CellField(data::Number, grid::Grid{T,AG}) where {T,AG}
  A = arraytype(grid)
  data = convert(A, fill(data, cell_size(grid)))
  return CellField(data, grid)
end

function FaceField(data::Number, grid::Grid{T,AG}) where {T,AG}
  A = arraytype(grid)
  data = convert(A, fill(data, face_size(grid)))
  return FaceField(data, grid)
end

FaceField(grid) = FaceField(0, grid)
CellField(grid) = CellField(0, grid)

CellField(data::Function, grid) = CellField(data.(grid.zc), grid)
FaceField(data::Function, grid) = FaceField(data.(grid.zf), grid)

#
# Basic 'Field' functionality
# 

@inline eachindex(c::CellField) = 1:c.grid.nz # interior indices of c
@inline eachindex(f::FaceField) = 2:f.grid.nz # interior indices of f

eltype(c::Field{A}) where A = eltype(A)

getindex(c::Field, inds...) = getindex(c.data, inds...)
setindex!(c::Field, d, inds...) = setindex!(c.data, d, inds...)
setindex!(c::Field, d::Field, inds...) = setindex!(c.data, d.data, inds...)

function setproperty!(sol::AbstractSolution, c::Symbol, data::Union{Number,AbstractArray,Function})
  set!(sol.c, data)
  return nothing
end

set!(c::Field{A}, data::Union{Number,AbstractArray}) where A = c.data .= convert(A, data)
set!(c::CellField{A}, data::Function) where A = c.data .= convert(A, data.(c.grid.zc))
set!(f::FaceField{A}, data::Function) where A = f.data .= convert(A, data.(f.grid.zf))
set!(c::Field{A,G}, d::Field{A,G}) where {A,G} = c.data .= convert(A, d.data)

similar(c::CellField) = CellField(c.grid)
similar(f::FaceField) = FaceField(f.grid)

# Define +, -, and * on fields as element-wise calculations on their data. This
# is only true for fields of the same type.
for TField in (:CellField, :FaceField)
    for op in (:+, :-, :*) 
        @eval begin
            # +, -, * a Field by a Number on the leTField.
            function $op(num::Number, f::$TField)
                ff = similar(f)
                @. ff.data = $op(num, f.data)
                ff  
            end 

            # +, -, * a Field by a Number on the right.
            $op(f::$TField, num::Number) = $op(num, f)

            # Multiplying two fields together
            function $op(f1::$TField, f2::$TField)
                f3 = similar(f1)
                @. f3.data = $op(f1.data, f2.data)
                f3  
            end 
        end 
    end 
end

#
# Differential operators defined on fields
# 

"Return the cell spacing at index i."
@inline dzc(c, i) = c.grid.dzc[i]
@inline dzc(c::Field{A,G}, i) where {A,G<:UniformGrid} = c.grid.dzc

"Return the face spacing at index i."
@inline dzf(c, i) = c.grid.dzf[i]
@inline dzf(c::Field{A,G}, i) where {A,G<:UniformGrid} = c.grid.dzf

@inline ∂z(a, i) = throw("∂z is not defined for arbitrary fields.")

"Return ∂c/∂z at index i."
@inline ∂z(c::CellField, i) = (c.data[i] - c.data[i-1]) / dzc(c, i)
@inline ∂²z(c::CellField, i) = (∂z(c, i+1) - ∂z(c, i)) / dzf(c, i)

"Return ∂f/∂z at index i."
@inline ∂z(f::FaceField, i) = (f.data[i+1] - f.data[i]) / dzf(f, i)
@inline ∂²z(f::FaceField, i) = (∂z(i, f) - ∂z(i-1, f)) / dzc(f)

∂z(a) = throw("∂z is not defined for arbitrary fields.")

"Calculate `f = ∂c/∂z` across the whole grid."
function ∂z!(f::FaceField, c::CellField)
  for i = eachindex(f)
    @inbounds f.data[i] = ∂z(c, i)
  end
  return nothing
end

"Calculate `c = ∂f/∂z` across the whole grid."
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
function ∂z(f::FaceField)
  c = CellField(f.grid)
  ∂z!(c, f)
  return c
end
