using
  Documenter,
  OceanTurb

makedocs(
     modules = [OceanTurb],
       clean = true,
     doctest = false,
   checkdocs = :all,
      format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
     authors = "Gregory L. Wagner",
    sitename = "OceanTurb.jl",

       pages = Any[
                "Home" => "index.md",
                "Numerical methods" => "numerics.md",
                "Boundary layer models" => Any[
                  "models/basics.md",
                  "models/kpp.md",
                  "models/pacanowskiphilander.md"],
                "DocStrings" => Any[
                      "man/types.md",
                      "man/functions.md"]
                      ]
)

deploydocs(
  repo = "github.com/glwagner/OceanTurb.jl.git",
)
