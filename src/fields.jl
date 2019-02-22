#= Fields for OceanTurb.jl

A `Field` is an abstractation that simplifies differential operations on a staggered grid.

The most import function that we output is `∂z`.
=#

import Base: setindex!, getindex, +, *, -

export
  Field,
  CellField,
  FaceField,
  dzc,
  dzf,
  ∂z,
  ∂z!

abstract type Field{A,G} end

struct CellField{A,G} <: Field{A,G} 
  data::A
  grid::G
end

struct FaceField{A,G} <: Field{A,G}
  data::A
  grid::G
end

CellField(grid::Grid{T,A}) where {T,A} = CellField(convert(A, zeros(cell_field_size(grid.nz))), grid)
FaceField(grid::Grid{T,A}) where {T,A} = FaceField(convert(A, zeros(face_field_size(grid.nz))), grid)

getindex(c::Field, inds...) = getindex(c.data, inds...)
setindex!(c::Field, d, inds...) = setindex!(c.data, d, inds...)
setindex!(c::Field, d::Field, inds...) = setindex!(c.data, d.data, inds...)

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
dzc(c, i) = c.grid.dzc[i]
dzc(c::Field{A,G}, i) where {A,G<:UniformGrid} = c.grid.dzc

"Return the face spacing at index i."
dzf(c, i) = c.grid.dzf[i]
dzf(c::Field{A,G}, i) where {A,G<:UniformGrid} = c.grid.dzf

∂z(a, i) = throw("∂z is not defined for arbitrary fields.")

"Return ∂c/∂z at index i."
∂z(c::CellField, i) = (c.data[i] - c.data[i-1]) / dzc(c, i)

"Return ∂f/∂z at index i."
∂z(f::FaceField, i) = (f.data[i+1] - f.data[i]) / dzf(c, i)

∂z(a) = throw("∂z is not defined for arbitrary fields.")

"Calculate `f = ∂c/∂z` across the whole grid."
function ∂z!(f::FaceField, c::CellField)
  for i = 1:c.grid.nz+1
    @inbounds f.data[i] = ∂z(c, i)
  end
  return nothing
end

"Calculate `c = ∂f/∂z` across the whole grid."
function ∂z!(c::CellField, f::FaceField)
  for i = 1:c.grid.nz
    @inbounds c.data[i] = ∂z(f, i)
  end
  return nothing
end

"Return the `FaceField` ∂c/∂z, where `c` is a `CellField`."
function ∂z(c::CellField{A,G}) where {A,G}
  f = FaceField(A, grid)
  ∂z!(f, c)
  return f
end

"Return the `CellField` ∂f/∂z, where `f` is a `FaceField`."
function ∂z(f::FaceField{A,G}) where {A,G}
  c = CellField(A, grid)
  ∂z!(c, f)
  return c
end
