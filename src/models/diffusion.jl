module Diffusion

using
    StaticArrays,
    OceanTurb

export
    Parameters,
    Model

import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top

# Just one field: "c"
@specify_solution CellField c

struct Parameters{T} <: AbstractParameters
    κ::T
end

struct Model{P, TS, G, E, T} <: AbstractModel{TS, G, E, T}
    @add_standard_model_fields
    parameters::P
end

function Model(; N=10, L=1.0, κ=0.1,
    grid = UniformGrid(N, L),
    parameters = Parameters(κ),
    stepper = :ForwardEuler,
    bcs = BoundaryConditions(ZeroFluxBoundaryConditions())
    )

    solution = Solution(CellField(grid))
    equation = Equation(calc_rhs_explicit!)
    timestepper = Timestepper(:ForwardEuler, solution)

    return Model(timestepper, grid, equation, solution, bcs, Clock(), parameters)
end

#
# Equation specification
#

# Equation specification
κ(m, i) = m.parameters.κ

function calc_rhs_explicit!(∂t, m)

    c = m.solution.c
    update_ghost_cells!(c, κ(m, 1), κ(m, c.grid.N), m, m.bcs.c)

    for i in eachindex(c)
        @inbounds ∂t.c[i] = ∇K∇c(κ(m, i+1), κ(m, i), c, i)
    end

    return nothing
end

end # module
