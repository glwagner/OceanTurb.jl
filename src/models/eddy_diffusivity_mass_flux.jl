module EDMF

using
    OceanTurb,
    StaticArrays,
    LinearAlgebra

import OceanTurb: ∇K∇c, ∇K∇c_bottom, ∇K∇c_top, Constants
import OceanTurb: oncell

import .KPP: ∂B∂z

const nsol = 5
@prefixed_solution ZeroPlume U V T S e

Base.@kwdef struct Parameters{T} <: AbstractParameters
         Cε :: T = 2.0    # Surface layer fraction
         Cκ :: T = 0.41   # Von Karman constant
         CK :: T = 0.1    # Minimum unresolved turbulence kinetic energy

    Ca_stab :: T = 2.7    # Stable buoyancy flux parameter for wind-driven turbulence
    Ca_unst :: T = -100.0 # Unstable buoyancy flux parameter for wind-driven turbulence

    Cn_stab :: T = 0.2    # Stable buoyancy flux parameter for wind-driven turbulence
    Cn_unst :: T = -1.0   # Unstable buoyancy flux parameter for wind-driven turbulence

        KU₀ :: T = 1e-6   # Interior viscosity for velocity
        KT₀ :: T = 1e-7   # Interior diffusivity for temperature
        KS₀ :: T = 1e-9   # Interior diffusivity for salinity
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

    KZP = ZeroPlumeAccessory{Function}(K, K, K, K, K)
    RZP = ZeroPlumeAccessory{Any}(RU, RV, nothing, nothing, Re)
    eqn = Equation(RZP, KZP, update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eqn, solution, lhs)

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

oncell(f::Function, m, i) = 0.5 * (f(m, i) + f(m, i+1))

turb_production(m, i) = K(m, i) * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2)
wb(m, i) = -K(m, i) * ∂B∂z(m, i)

#
# Equation entry
#
K(m, i) = m.parameters.CK * mixing_length(m, m.grid.zf[i], i) * sqrt(m.solution.e[i])

De(m, i) = m.parameters.Cε * m.solution.e[i]^(3/2) / mixing_length(m, m.grid.zc[i], i)

RU(m, i) =   m.constants.f * m.solution.V[i]
RV(m, i) = - m.constants.f * m.solution.U[i]
Re(m, i) = oncell(turb_production, m, i) + oncell(wb, m, i) - De(m, i)

end # module
