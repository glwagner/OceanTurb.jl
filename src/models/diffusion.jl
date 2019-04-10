module Diffusion

using LinearAlgebra, OceanTurb
import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top

# Just one field: "c"
@solution c

struct Parameters{T} <: AbstractParameters
    K::T
end

struct Model{P, TS, G, T} <: AbstractModel{TS, G, T}
    @add_standard_model_fields
    parameters::P
end

function Model(; N=10, L=1.0, K=0.1,
          grid = UniformGrid(N, L),
    parameters = Parameters(K),
       stepper = :ForwardEuler,
           bcs = BoundaryConditions(ZeroFluxBoundaryConditions())
    )

    solution = Solution(CellField(grid))
    K = Accessory(Kc)
    R = Accessory(nothing)
    eqn = Equation(R, K)
    lhs = LeftHandSide(solution)
    timestepper = Timestepper(stepper, eqn, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, parameters)
end

#
# Equation specification
#

Rc(m, i) = nothing
Kc(m, i) = m.parameters.K

end # module
