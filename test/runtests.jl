using
    OceanTurb,
    LinearAlgebra,
    Test

Timer(t -> println(" "), 0, interval=240)

#
# Run tests
#

steppers = (:ForwardEuler, :BackwardEuler)

include("test_utils.jl")
include("test_solvers.jl")
include("test_grids.jl")
include("test_fields.jl")

# Models
include("test_diffusion.jl")
include("test_pp.jl")
include("test_kpp.jl")
include("test_modularkpp.jl")
include("test_tkemassflux.jl")
