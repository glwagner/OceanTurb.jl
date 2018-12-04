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

abstract type AbstractParameters end
abstract type ModelParameters end
abstract type Model end

# Gregorian calendar
const second = 1.0
const minute = 60second
const hour = 60minute
const day = 24hour
const year = 365day
# Rotation rate
const stellaryear = 23hour + 56minute + 4.098903691
const Ω = 2π/stellaryear

include("utils.jl")
include("forcing.jl")
include("parameters.jl")
include("ocean.jl")

include("PriceWellerPinkel/PriceWellerPinkel.jl")

end # module
