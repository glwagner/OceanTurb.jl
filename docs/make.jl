using
  Documenter,
  OceanTurb

# Set up a timer to print a space ' ' every 240 seconds. This is to avoid Travis CI
# timing out when building demanding Literate.jl examples.
Timer(t -> println(" "), 0, interval=240)

makedocs(
     modules = [OceanTurb],
       clean = true,
     doctest = false,
   checkdocs = :all,
      format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                               mathengine = Documenter.MathJax2()),
     authors = "Gregory L. Wagner",
    sitename = "OceanTurb.jl",

       pages = Any[
                "Home" => "index.md",
                "Turbulence, fluxes, and physics" =>  "basics.md",
                "Numerical methods" => "numerics.md",
                "Turbulence models" => Any[
                  "models/kpp.md",
                  "models/modular_kpp.md",
                  "models/tke_mass_flux.md",
                  "models/pacanowskiphilander.md"],
                "DocStrings" => Any[
                      "man/types.md",
                      "man/functions.md"]
                      ]
)

deploydocs(        repo = "github.com/glwagner/OceanTurb.jl.git",
              versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
          push_preview = true
          )
