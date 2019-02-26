#= Fields for OceanTurb.jl

A `Field` is an abstractation that simplifies differential operations on a staggered grid.

OceanTurb.jl solves one-dimensional PDEs on a staggered grid with ghost cells.
The geometry of a grid with nz=3 is 

                                        ^ z 
                                        |   
        
                      ---
         i=4           *    
                j=4   ===       
         i=3           *    
                j=3   ---
         i=2           *    ^            
                j=2   ---   | dzc (j=2) 
         i=1           *    v  
                j=1   ===     ^        
         i=0           *      | dzf  (i=0)
                      ---     v
              

where the i's index cells and the j's index faces. The variable dzc gives the separatiol between
cell centers, and dzf gives the separation between faces.

Accordingly, variables located in cells (`CellFields`) have dimension nz+2, 
and variables located at cell faces (`FaceFields`) have dimension dimension nz+1.

The most import function that we output is `∂z`.
=#

import Base: setindex!, getindex, eachindex, +, *, -, similar

export
  Field,
  CellField,
  FaceField,
  dzc,
  dzf,
  ∂z,
  ∂²z,
  ∂z!

abstract type Field{T,A,G} end

struct CellField{T,A,G} <: Field{T,A,G} 
  data::OffsetArray{T,1,A}
  grid::G
end

struct FaceField{T,A,G} <: Field{T,A,G}
  data::A
  grid::G
end

function CellField(data::Array, grid::Grid{T,A}) where {T,A}
  sz = cell_field_size(grid.nz)
  data = OffsetVector(A(data), 0:sz-1)
  return CellField(data, grid)
end

function CellField(grid::Grid{T,A}) where {T,A} 
  sz = cell_field_size(grid.nz)
  data = OffsetVector(zeros(eltype(A), sz), 0:sz-1)
  return CellField(data, grid)
end

function FaceField(grid::Grid{T,A}) where {T,A} 
  sz = face_field_size(grid.nz)
  data = zeros(eltype(A), sz)
  return FaceField{T,A,typeof(grid)}(data, grid)
end

FaceField(data::Array, grid::Grid{T,A}) where {T,A} = FaceField{T,A,typeof(grid)}(data, grid)

@inline eachindex(c::CellField) = 1:c.grid.nz
@inline eachindex(f::FaceField) = 1:f.grid.nz+1

getindex(c::Field, inds...) = getindex(c.data, inds...)
setindex!(c::Field, d, inds...) = setindex!(c.data, d, inds...)
setindex!(c::Field, d::Field, inds...) = setindex!(c.data, d.data, inds...)

similar(c::CellField) = CellField(c.data, c.grid)
similar(f::FaceField) = FaceField(f.data, f.grid)

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
@inline dzc(c::Field{T,A,G}, i) where {T,A,G<:UniformGrid} = c.grid.dzc

"Return the face spacing at index i."
@inline dzf(c, i) = c.grid.dzf[i]
@inline dzf(c::Field{T,A,G}, i) where {T,A,G<:UniformGrid} = c.grid.dzf

@inline ∂z(a, i) = throw("∂z is not defined for arbitrary fields.")

"Return ∂c/∂z at index i."
@inline ∂z(c::CellField, i) = (c.data[i] - c.data[i-1]) / dzc(c, i)
@inline ∂²z(c::CellField, i) = (∂z(i+1, c) - ∂z(i, c)) / dzf(c)

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
function ∂z(c::CellField{T,A,G}) where {T,A,G}
  f = FaceField(c.grid)
  ∂z!(f, c)
  return f
end

"Return the `CellField` ∂f/∂z, where `f` is a `FaceField`."
function ∂z(f::FaceField{T,A,G}) where {T,A,G}
  c = CellField(f.grid)
  ∂z!(c, f)
  return c
end
