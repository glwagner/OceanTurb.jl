"Example model forced by Southern Ocean data."
function loadexample(H=400, nz=400)
  filename = "example_southern_ocean_forcing.jld2"
  datapath = joinpath(dirname(pathof(OceanMixedLayerModels)), "..", "data")
  filepath = joinpath(datapath, filename)
  forcing = loadforcing(filepath)
  PriceWellerPinkelModel(forcing=forcing, ocean=Ocean(H=H, nz=nz))
end
