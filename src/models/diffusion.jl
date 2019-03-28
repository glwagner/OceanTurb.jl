module Diffusion

using
    StaticArrays,
    LinearAlgebra,
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

struct Model{P, TS, G, T} <: AbstractModel{TS, G, T}
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

    if implicit(stepper)
        get_κ = SolutionLike(κ)
        lhs_array = OceanTurb.build_lhs_array(solution)
        lhs = SolutionLike(lhs_array)
        timestepper = Timestepper(stepper, calc_rhs_implicit!, solution, lhs)
    else
        timestepper = Timestepper(stepper, calc_rhs_explicit!, solution)
    end

    return Model(Clock(), grid, timestepper, solution, bcs, parameters)
end

#
# Equation specification
#

κ(m, i) = m.parameters.κ

function calc_rhs_explicit!(∂t, m)
    c = m.solution.c
    update_ghost_cells!(c, κ(m, 1), κ(m, c.grid.N), m, m.bcs.c)
    for i in eachindex(c)
        @inbounds ∂t.c[i] = ∇K∇c(κ(m, i+1), κ(m, i), c, i)
    end

    return nothing
end

function calc_rhs_implicit!(rhs, m)
    c = m.solution.c
    update_ghost_cells!(c, κ(m, 1), κ(m, c.grid.N), m, m.bcs.c)

    # Add flux across top and bottom boundary
    @inbounds begin
        rhs.c[m.grid.N] = top_flux_div(0, κ(m, m.grid.N), c)
        rhs.c[1] = bottom_flux_div(0, κ(m, 1), c)
    end
    return nothing
end

end # module
