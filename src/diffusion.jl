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
    equation = Equation(calc_rhs!)
    timestepper = Timestepper(:ForwardEuler, solution)

    return Model(timestepper, grid, equation, solution, bcs, Clock(), parameters)
end

#
# Equation specification
#

# Equation specification
κ(m, i) = m.parameters.κ

function calc_rhs!(∂t, m)

    c = m.solution.c

    @inbounds begin
        for i in interior(c)
            @inbounds ∂t.c[i] = ∇K∇c(κ(m, i+1), κ(m, i), c, i)
        end

        i = m.grid.N
        ∂t.c[i] = ∇K∇c_top(κ(m, i+1), κ(m, i), c, m.bcs.c.top, m)

        i = 1
        ∂t.c[i] = ∇K∇c_bottom(κ(m, i+1), κ(m, i), c, m.bcs.c.bottom, m)
    end

    return nothing
end

end # module
