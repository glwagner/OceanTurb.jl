module OceanTurb

export
  Model,
  Grid,
  RegularGrid
    
using
    JLD2,
    Interpolations

using Statistics: mean

# Gregorian calendar-ic globals
const second = 1.0
const minute = 60second
const hour = 60minute
const day = 24hour
const year = 365day
const stellaryear = 23hour + 56minute + 4.098903691
const Ω = 2π/stellaryear

include("utils.jl")

abstract type Grid{T} end
abstract type Variables{T} end
abstract type Parameters{T} end

mutable struct Model{T,A}
  grid::Grid{T}
  variables::Variables{A}
  parameters::Parameters{T}
  t::T
  step::Int
end

Model(nz, H, κ) = Model{Float64,Array{Float64,1}}(
  RegularGrid(nz, H), PrimitiveVariables(nz), DiffusionParameters(κ), 0.0, 0)

struct DiffusionParameters{T} <: Parameters{T}
  κ::T
end

struct PrimitiveVariables{A} <: Variables{A}
  U::A
  V::A
  T::A
  S::A
end

PrimitiveVariables(nz) = PrimitiveVariables(zeros(nz), zeros(nz), zeros(nz), zeros(nz))

"""
    RegularGrid(T, H, nz)

Construct a 1D grid.
"""
struct RegularGrid{T,A} <: Grid{T}
  nz::Int
  dz::T
  zC::A
  zE::A
end

function RegularGrid(T, nz, H)
  dz = H/nz
  zE = collect(T, range(-H, step=dz, stop=0)) # cell faces
  zC = 0.5*(zE[1:end-1] + zE[2:end]) # cell centers
  RegularGrid(nz, dz, zC, zE)
end

RegularGrid(nz, H) = RegularGrid(Float64, nz, H)


end # module
