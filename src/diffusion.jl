module Diffusion

using
    StaticArrays,
    OceanTurb

export
    Parameters,
    Model

# Just one field: "c"
@specify_solution CellField c

struct Parameters{T} <: AbstractParameters
    κ::T
end

struct Model{PT, TS, G, E, T} <: AbstractModel{TS, G, E, T}
    @add_standard_model_fields
    parameters::Parameters{PT}
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

# Convenient operators
κ∂z(κ::Number, c, i) = κ*∂z(c, i)
κ∂z(κ::AbstractField, c, i) = κ.data[i]*∂z(c, i)
κ∂z(κ::Function, c, i) = κ(c.grid.zf[i]) * ∂z(c, i) # works for κ(z)

∇κ∇c(κ, c, i)           = ( κ∂z(κ, c, i+1) -    κ∂z(κ, c, i)      ) /    Δf(c, i)
∇κ∇c_top(κ, c, flux)    = (     -flux      - κ∂z(κ, c, length(c)) ) / Δf(c, length(c))
∇κ∇c_bottom(κ, c, flux) = (  κ∂z(κ, c, 2)  +        flux          ) /    Δf(c, 1)

# Boundary conditions...
const BC = BoundaryCondition

# Top and bottom flux estimates for constant (Dirichlet) boundary conditions
bottom_flux(κ, c, c_bndry, Δf) = -2*bottom(κ)*( bottom(c) - c_bndry ) / bottom(Δf) # -κ*∂c/∂z at the bottom
top_flux(κ, c, c_bndry, Δf)    = -2*  top(κ) *(  c_bndry  -  top(c) ) /   top(Δf)  # -κ*∂c/∂z at the top

∇κ∇c_bottom(m, bc::BC{<:Flux}) = ∇κ∇c_bottom(m.parameters.κ, m.solution.c, getbc(m, bc))
∇κ∇c_top(m, bc::BC{<:Flux}) = ∇κ∇c_top(m.parameters.κ, m.solution.c, getbc(m, bc))

function ∇κ∇c_bottom(model, bc::BC{<:Value})
    flux = bottom_flux(model.parameters.κ, model.solution.c, getbc(model, bc), model.grid.Δf)
    return ∇κ∇c_bottom(model.parameters.κ, model.solution.c, flux)
end

function ∇κ∇c_top(model, bc::BC{<:Value})
    flux = top_flux(model.parameters.κ, model.solution.c, getbc(model, bc), model.grid.Δf)
    return ∇κ∇c_top(model.parameters.κ, model.solution.c, flux)
end

function calc_rhs!(rhs, model)
    for i in interior(model.solution.c)
        @inbounds rhs.c.data[i] = ∇κ∇c(model.parameters.κ, model.solution.c, i)
    end

    rhs.c.data[end] = ∇κ∇c_top(model, model.bcs.c.top)
    rhs.c.data[1] = ∇κ∇c_bottom(model, model.bcs.c.bottom)

    return nothing
end

end # module
