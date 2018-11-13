module PriceWellerPinkel

export
  Forcing,
  Profile,
  Model,
  density,
  convect!,
  @zeros,
  loadexample,
  TestModel

using
  JLD2

using Statistics: mean

const DEBUG = true

const g = 9.807
const Cₚ = 3900.0

macro zeros(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end

include("forcing.jl")
include("model.jl")
include("eos.jl")
include("mixing.jl")

function loadexample(H=400, nz=400)
  datapath = joinpath(dirname(pathof(PriceWellerPinkel)), "..", "data")
  filename = "example_southern_ocean_forcing.jld2"
  filepath = joinpath(datapath, filename)
  forcing = loadforcing(filepath)
  Model(forcing, H, nz) 
end

TestModel(H, nz, nt=4) = Model(Forcing(zeros(nt)), H, nz) 

end # module
