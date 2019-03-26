#= Grid types for OceanTurb.jl.

OceanTurb.jl solves one-dimensional PDEs on a staggered grid.
The geometry of a grid with N=3 is

```
      ▲ z
      |

         i=4           *    (ghost cell)
                j=4   ===       ▲
         i=3           *        | Δf[3]
                j=3   ---       ▼
         i=2           *    ▲
                j=2   ---   | Δc[2]
         i=1           *    ▼
                j=1   ===
         i=0           *    (ghost cell)
```

where the i's index cells and the j's index faces.
The variable Δc gives the separation between
cell centers, and Δf gives the separation between faces.

Accordingly, variables located in cells (`CellFields`) have dimension N,
and variables located at cell faces (`FaceFields`) have dimension dimension N+1.
=#

import Base: eltype, length, size

# Staggered grid lengths and sizes for fields
cell_length(N) = N+2
face_length(N) = N+1
  cell_size(N) = (N+2,)
  face_size(N) = (N+1,)

     height(g::Grid) = g.L
     length(g::Grid) = g.N
       size(g::Grid) = (g.N,)
cell_length(g::Grid) = cell_length(g.N)
face_length(g::Grid) = face_length(g.N)
  cell_size(g::Grid) = cell_size(g.N)
  face_size(g::Grid) = face_size(g.N)

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
    UniformGrid([A, T], L, N)

Construct a 1D finite-volume grid with type `T` and array type `A`,
with `N` cells on the domain `z = [0, L]`.
A `Grid` has two type parameters: an element type `T`,
and and an array type `A`. `T` and `A` default to `Float64` and `Array{T,1}`,
respectively.
"""
struct UniformGrid{T, A} <: Grid{T, A}
  N  :: Int
  L  :: T
  Δc :: T
  Δf :: T
  zc :: A
  zf :: A
end

function UniformGrid(T, N, L)
  Δ = convert(T, L/N)
  L = convert(T, L)
  half_Δ = convert(T, 0.5Δ)

  zc = range(-L+half_Δ; length=N, stop=-half_Δ)
  zf = range(-L; length=N+1, stop=zero(T))

  UniformGrid(N, L, Δ, Δ, zc, zf)
end

# Defaults
UniformGrid(N=3, L=1) = UniformGrid(Float64, N, L)

"Return the cell spacing at index i."
Δc(grid::UniformGrid, i) = grid.Δc

"Return the face spacing at index i."
Δf(grid::UniformGrid, i) = grid.Δf
