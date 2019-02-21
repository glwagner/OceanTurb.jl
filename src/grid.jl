"""
    Grid(T, H, nz)

Construct a 1D grid.
"""
struct Grid{T}
  nz::Int
  dz::T
  zC::AbstractArray{T}
  zE::AbstractArray{T}
end

function Grid(; Tel=Float64, H=400, nz=100)
  dz = H/nz
  zE = collect(Tel, range(-H, step=dz, stop=0)) # cell faces
  zC = 0.5*(zF[1:end-1] + zF[2:end]) # cell centers
  Grid(nz, dz, zC, zE)
end

Grid(H, nz) = Grid(Float64, H, nz)
