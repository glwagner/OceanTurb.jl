module OceanBoundaryLayerModels

export
    Ocean,

    AbstractParameters,
    ModelParameters,
    Parameters,

    Forcing,
    ForcingInterpolant,

    Model,

    density,
    stepforward!,
    updatevars!,
    mixedlayerdepth,

    loadforcing,

    Ω,

    @zeros

using
    JLD2,
    Interpolations

using Statistics: mean

const DEBUG = true

abstract type Model end
abstract type AbstractParameters end

# Gregorian calendar-ic globals
const second = 1.0
const minute = 60second
const hour = 60minute
const day = 24hour
const year = 365day
const stellaryear = 23hour + 56minute + 4.098903691
const Ω = 2π/stellaryear

include("utils.jl")
include("forcing.jl")

struct Constants{T}
  f::T
  g::T
  ρ0::T
  ν::T
  κT::T
  κS::T
end

struct Grid{T,A}
  nz::Int
  dz::T
  zC::A
  zF::A
end

struct State{A}
  u::A
  v::A
  T::A
  S::A
  ρ::A
end

mutable struct Clock{T}
  step::Int
  t::Float64
end

struct Problem # ?
  consts::Constants
  forcing::Forcing
  grid::Grid
  state::State
  model::Model
    #vars
    #params
  clock::Clock
  timestepper
end

struct PhilanderPacanowski <: Model
end


include("idealizedforcing.jl")

#include("parameters.jl")
#include("ocean.jl")
#include("PriceWellerPinkel/PriceWellerPinkel.jl")

end # module
