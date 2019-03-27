module SouzaAdjustment

using
    StaticArrays,
    OceanTurb

export
    Parameters,
    Model

import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top

@pair_specify_solution CellField B FaceField wb

struct Parameters{T} <: AbstractParameters
    Cκ :: T  # Diffusivity
    Cs :: T  # Stability - flux proportionality
    Cc :: T  # Cut-off in stability - flux relationship
    Cτ :: T  # Relaxation time-scale for buoyancy flux
end

function Parameters(T=Float64;
     κ = 1e-5,
    Cs = 1.0,
    Cc = 0.0,
    Cτ = 3600
    )

    Parameters{T}(κ, Cs, Cc, Cτ)
end

struct Model{TS, G, E, T} <: AbstractModel{TS, G, E, T}
    @add_standard_model_fields
    parameters::Parameters
end

function Model(; N=10, L=1.0,
    grid = UniformGrid(N, L),
    parameters = Parameters(),
    stepper = :ForwardEuler,
    bcs = BoundaryConditions(ZeroFluxBoundaryConditions(), ZeroFluxBoundaryConditions())
    )

    solution = Solution(CellField(grid), FaceField(grid))
    equation = Equation(calc_rhs_explicit!)
    timestepper = Timestepper(:ForwardEuler, solution)

    return Model(timestepper, grid, equation, solution, bcs, Clock(), parameters)
end

#
# Equation specification
#

# Equation specification
κ(m, i) = m.parameters.Cκ

function calc_rhs_explicit!(∂t, m)

    B, wb = m.solution
    κ, Cs, Cc, τ = m.parameters.Cκ, m.parameters.Cs, m.parameters.Cc, m.parameters.Cτ

    update_ghost_cells!(B, κ, κ, m, m.bcs.c)

    for i in eachindex(B)
        @inbounds ∂t.B[i] = -∂z(wb, i) + κ*∂²z(B, i)
    end

    for i in interiorindices(wb)
        @inbounds ∂t.wb[i] = κ*∂²z(wb, i) - Cs*min(Cc, ∂z(B, i)) - 1/τ * wb[i]
    end

    return nothing
end

end # module
