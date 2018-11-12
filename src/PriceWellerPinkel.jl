module PriceWellerPinkel

export
  Forcing,
  Profile,
  Model,
  @zeros,
  loadexample

using
  JLD2

macro zeros(T, dims, vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
  expr
end

include("forcing.jl")
include("model.jl")

function loadexample(H=400, nz=400)
  datapath = joinpath(dirname(pathof(PriceWellerPinkel)), "..", "data")
  filename = "example_southern_ocean_forcing.jld2"
  filepath = joinpath(datapath, filename)
  forcing = loadforcing(filepath)
  Model(forcing, H, nz) 
end


end # module
