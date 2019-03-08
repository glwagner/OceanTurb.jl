#= Grid types for OceanTurb.jl.

OceanTurb.jl solves one-dimensional PDEs on a staggered grid.
The geometry of a grid with nz=3 is

```
      ▲ z
      |

                j=4   ===       ▲
         i=3           *        | dzf[3]
                j=3   ---       ▼
         i=2           *    ▲
                j=2   ---   | dzc[2]
         i=1           *    ▼
                j=1   ===
```

where the i's index cells and the j's index faces.
The variable dzc gives the separation between
cell centers, and dzf gives the separation between faces.

Accordingly, variables located in cells (`CellFields`) have dimension nz,
and variables located at cell faces (`FaceFields`) have dimension dimension nz+1.
=#

import Base: eltype, length, size

# Staggered grid lengths and sizes for fields
cell_length(nz) = nz
face_length(nz) = nz+1
  cell_size(nz) = (nz,)
  face_size(nz) = (nz+1,)

     height(g::Grid) = g.Lz
     length(g::Grid) = g.nz
       size(g::Grid) = (g.nz,)
cell_length(g::Grid) = cell_length(g.nz)
face_length(g::Grid) = face_length(g.nz)
  cell_size(g::Grid) = cell_size(g.nz)
  face_size(g::Grid) = face_size(g.nz)

height(m::AbstractModel) = height(m.grid)
length(m::AbstractModel) = length(m.grid)
  size(m::AbstractModel) = size(m.grid)

eltype(::Grid{T}) where T = T

"""
    arraytype(grid::Grid)

Return the array type corresponding to data
that lives on `grid`. Defaults to `Array`.
New data types (for example, grids that exist on GPUs) must
implement new array types.
"""
arraytype(grid::Grid{T}) where T = Array{T,1} # default array type

"""
    UniformGrid([A, T], Lz, nz)

Construct a 1D finite-volume grid with type `T` and array type `A`,
with `nz` cells on the domain `z = [0, Lz]`.
A `Grid` has two type parameters: an element type `T`,
and and an array type `A`. `T` and `A` default to `Float64` and `Array{T,1}`,
respectively.
"""
struct UniformGrid{T,A} <: Grid{T,A}
  nz::Int
  Lz::T
  dzc::T
  dzf::T
  zc::A
  zf::A
end

function UniformGrid(T, nz, Lz)
  dz = convert(T, Lz/nz)
  Lz = convert(T, Lz)
  half_dz = convert(T, 0.5dz)

  zc = range(-Lz+half_dz; length=cell_length(nz), stop=-half_dz)
  zf = range(-Lz; length=face_length(nz), stop=zero(T))

  UniformGrid(nz, Lz, dz, dz, zc, zf)
end

# Defaults
UniformGrid(nz=1, Lz=1.0) = UniformGrid(Float64, nz, Lz)
