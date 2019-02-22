#= Grid types for OceanTurb.jl.

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
=#

export
  Grid,
  UniformGrid

abstract type Grid{T,A} end

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

function UniformGrid(A, T, nz, Lz)
  dz = Lz/nz
  zc = collect(T, range(-Lz-0.5*dz, step=dz, stop=0.5*dz)) # centers with ghost cells
  zf = collect(T, range(-Lz, step=dz, stop=0)) # faces
  UniformGrid{T,A}(nz, Lz, dz, dz, zc, zf)
end

# Defaults
UniformGrid(T, nz, Lz) = UniformGrid(Array{T,1}, T, nz, Lz)
UniformGrid(nz, Lz) = UniformGrid(Float64, nz, Lz)

# Staggered grid sizes for fields
cell_field_size(nz) = nz+2
face_field_size(nz) = nz+1
