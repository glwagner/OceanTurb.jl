module PriceWellerPinkel

export
  Forcing,
  ForcingInterpolant,
  Profile,
  Model,
  density,
  convect!,
  @zeros,
  loadexample

using
  JLD2,
  Interpolations

using Statistics: mean

const DEBUG = true
const g = 9.807

# Gregorian calendar
const second = 1.0
const minute = 60second
const hour = 60minute
const day = 24hour
const year = 365day

# Rotation rate
const stellaryear = 23hour + 56minute + 4.098903691
const Ω = 2π/stellaryear

macro zeros(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end

include("forcing.jl")
include("model.jl")
include("eos.jl")
include("mixing.jl")

"Example model forced by Southern Ocean data."
function loadexample(H=400, nz=400)
  datapath = joinpath(dirname(pathof(PriceWellerPinkel)), "..", "data")
  filename = "example_southern_ocean_forcing.jld2"
  filepath = joinpath(datapath, filename)
  forcing = loadforcing(filepath)
  Model(forcing=forcing, H=H, nz=nz) 
end

end # module
