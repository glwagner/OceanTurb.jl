module Diffusion

using OceanTurb

# Just one field: "c"
@solution c

Base.@kwdef struct Parameters{T} <: AbstractParameters
    K :: T = 0.0
    W :: T = 0.0
    μ :: T = 0.0
end

struct Model{P, TS, G, T} <: AbstractModel{TS, G, T}
    clock       :: Clock{T}
    grid        :: G
    timestepper :: TS
    solution    :: Solution
    bcs         :: BoundaryConditions
    parameters  :: P
end

function Model(; N=10, L=1.0, K=0.1, W=0.0, μ=0.0,
          grid = UniformGrid(N, L),
    parameters = Parameters(K, W, μ),
       stepper = :ForwardEuler,
           bcs = BoundaryConditions(ZeroFluxBoundaryConditions())
    )

    solution = Solution(CellField(grid))

      K = (Kc,)
      W = (Wc,)
      L = (Lc,)

    eqn = Equation(K=K, M=W, L=L)
    lhs = build_lhs(solution)

    timestepper = Timestepper(stepper, eqn, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, parameters)
end

@inline Kc(m, i) = m.parameters.K
@inline Wc(m, i) = m.parameters.W
@inline Lc(m, i) = m.parameters.μ

end # module
