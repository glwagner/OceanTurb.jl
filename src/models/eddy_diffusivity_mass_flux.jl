module EDMF

using
    OceanTurb,
    StaticArrays,
    LinearAlgebra

import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top, Constants
import OceanTurb: oncell

import .KPP: ∂B∂z

const nsol = 5
@named_solution ZeroPlume U V T S e

struct Parameters{T} <: AbstractParameters
         Cε :: T  # Surface layer fraction
         Cκ :: T  # Von Karman constant
         CK :: T  # Minimum unresolved turbulence kinetic energy

    Ca_stab :: T  # Stable buoyancy flux parameter for wind-driven turbulence
    Ca_unst :: T  # Unstable buoyancy flux parameter for wind-driven turbulence

    Cn_stab :: T  # Stable buoyancy flux parameter for wind-driven turbulence
    Cn_unst :: T  # Unstable buoyancy flux parameter for wind-driven turbulence

        KU₀ :: T  # Interior viscosity for velocity
        KT₀ :: T  # Interior diffusivity for temperature
        KS₀ :: T  # Interior diffusivity for salinity
end

function Parameters(T=Float64;
         Cε = 2.0,
         Cκ = 0.41,
         CK = 0.1,
    Ca_stab = 2.7,
    Ca_unst = -100,
    Cn_stab = 0.2,
    Cn_unst = -1,
        KU₀ = 1e-6,
        KT₀ = 1e-7,
        KS₀ = 1e-9
    )

    Parameters{T}(Cε, Cκ, CK, Ca_stab, Ca_unst, Cn_stab, Cn_unst, KU₀, KT₀, KS₀)
end

mutable struct State{T} <: FieldVector{5, T}
    Fu :: T
    Fv :: T
    Fθ :: T
    Fs :: T
    Fb :: T
end

State(T=Float64) = State{T}(0, 0, 0, 0, 0)

"""
    update_state!(model)

Update the top flux conditions and mixing depth for `model`
and store in `model.state`.
"""
function update_state!(m)
    m.state.Fu = getbc(m, m.bcs.U.top)
    m.state.Fv = getbc(m, m.bcs.V.top)
    m.state.Fθ = getbc(m, m.bcs.T.top)
    m.state.Fs = getbc(m, m.bcs.S.top)
    m.state.Fb = m.constants.g * (m.constants.α * m.state.Fθ - m.constants.β * m.state.Fs)
    return nothing
end

struct ZeroPlumeModel{TS, G, T} <: AbstractModel{TS, G, T}
    @add_clock_grid_timestepper
      solution :: ZeroPlumeSolution
           bcs :: ZeroPlumeBoundaryConditions
    parameters :: Parameters{T}
    constants  :: Constants{T}
         state :: State{T}
end

function ZeroPlumeModel(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
      parameters = Parameters(),
         stepper = :ForwardEuler,
             bcs = ZeroPlumeBoundaryConditions((ZeroFluxBoundaryConditions() for i=1:nsol)...)
    )

    solution = ZeroPlumeSolution((CellField(grid) for i=1:nsol)...)
    timestepper = Timestepper(stepper, calc_rhs_explicit!, solution)

    return ZeroPlumeModel(Clock(), grid, timestepper, solution, bcs, parameters, constants, State())
end


# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

mixing_length(Cκ, Ca, Cn, Fb, ωτ, z) = Cκ * z * (1 - z * Ca * Fb / ωτ^3)^Cn

function mixing_length(m, z, i)
    if isunstable(m)
        mixing_length(m.parameters.Cκ,
                      m.parameters.Ca_unst, m.parameters.Cn_unst, ωτ(m), z)
    else
        mixing_length(m.parameters.Cκ,
                      m.parameters.Ca_stab, m.parameters.Cn_stab, ωτ(m), z)
    end
    return nothing
end

K(m, i) = m.parameters.CK * mixing_length(m, m.grid.zf[i], i) * sqrt(m.solution.e[i])

oncell(f::Function, m, i) = 0.5 * (f(m, i) + f(m, i+1))

turb_production(m, i) = K(m, i) * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2)
wb(m, i) = -K(m, i) * ∂B∂z(m, i)

#
# Equation entry
#

function calc_rhs_explicit!(∂t, m)

    # Preliminaries
    update_state!(m)
    U, V, T, S, e = m.solution

    N = m.grid.N
    update_ghost_cells!(U, K(m, 1), K(m, N), m, m.bcs.U)
    update_ghost_cells!(V, K(m, 1), K(m, N), m, m.bcs.V)
    update_ghost_cells!(T, K(m, 1), K(m, N), m, m.bcs.T)
    update_ghost_cells!(S, K(m, 1), K(m, N), m, m.bcs.S)
    update_ghost_cells!(e, K(m, 1), K(m, N), m, m.bcs.e)

    for i in eachindex(U)
        @inbounds begin
            ∂t.U[i] = ∇K∇c(K(m, i+1), K(m, i), U, i) + m.constants.f * V[i]
            ∂t.V[i] = ∇K∇c(K(m, i+1), K(m, i), V, i) - m.constants.f * U[i]
            ∂t.T[i] = ∇K∇c(K(m, i+1), K(m, i), T, i)
            ∂t.S[i] = ∇K∇c(K(m, i+1), K(m, i), S, i)
            ∂t.e[i] = (
                ∇K∇c(K(m, i+1), K(m, i), e, i)
                + oncell(turb_production, m, i) + oncell(wb, m, i)
                - m.parameters.Cε * e[i]^(3/2) / mixing_length(m, m.grid.zc[i], i))
        end
    end

    return nothing
end

end # module
