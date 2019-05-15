module Diffusion

using LinearAlgebra, OceanTurb
import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top
import Base: @propagate_inbounds

# Just one field: "c"
@solution c

struct Parameters{T} <: AbstractParameters
    K::T
end

struct Model{P, TS, G, T} <: AbstractModel{TS, G, T}
    clock       :: Clock{T}
    grid        :: G
    timestepper :: TS
    solution    :: Solution
    bcs         :: BoundaryConditions
    parameters  :: P
end

function Model(; N=10, L=1.0, K=0.1,
          grid = UniformGrid(N, L),
    parameters = Parameters(K),
       stepper = :ForwardEuler,
           bcs = BoundaryConditions(ZeroFluxBoundaryConditions())
    )

    solution = Solution(CellField(grid))
    K = Accessory(Kc)
    R = Accessory(Rc) #nothing)
    eqn = Equation(R, K)
    lhs = build_lhs(solution) #LeftHandSide(solution)
    timestepper = Timestepper(stepper, eqn, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, parameters)
end

#
# Equation specification
#

@inline Rc(m, i) = 0.0
@inline Kc(m, i) = m.parameters.K

end # module
