using Pkg
Pkg.activate("..")

using Documenter, OceanMixedLayerModels

makedocs(
   modules = [OceanMixedLayerModels],
   doctest = false, 
     clean = true,
 checkdocs = :all,
    format = :html,
   authors = "Gregory L. Wagner",
  sitename = "OceanMixedLayerModels.jl",
     pages = Any[
              "Home" => "index.md",
                "DocStrings" => Any[
                    "man/types.md",
                    "man/functions.md"]
                    ]
)

#=
deploydocs(
       repo = "github.com/glwagner/PriceWellerPinkel.jl.git",
     target = "build",
      julia = "1.0",
     osname = "linux",
       deps = nothing,
       make = nothing
 )
 =#
