#= Grid types for OceanTurb.jl.

See fields.jl for the geometry of an OceanTurb.Field.
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
arraytype(grid::Grid{T}) where T = Array{T, 1} # default array type

"""
    UniformGrid([T], L, N)

Construct a 1D finite-volume grid with type `T` and array type `A`,
with `N` cells on the domain `z = [0, L]`.
A `Grid` has two type parameters: an element type `T`,
and and a range type `R`. `T` defaults to Float64.
"""
struct UniformGrid{T, R} <: Grid{T, R}
    N  :: Int
    L  :: T
    Δc :: T
    Δf :: T
    zc :: R
    zf :: R
end

function UniformGrid(T, N::Int, L::Number)
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
@inline Δc(grid::UniformGrid, i) = grid.Δc

"Return the face spacing at index i."
@inline Δf(grid::UniformGrid, i) = grid.Δf
