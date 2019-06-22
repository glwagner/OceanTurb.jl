module KPP_TKE

using OceanTurb

import ..OceanTurb: oncell, onface
import .KPP: ∂B∂z, isunstable, ωτ, ωb
import .ModularKPP: AbstractModularKPPModel

const nsol = 5
@solution U V T S e

Base.@kwdef struct TKEParameters{T} <: AbstractParameters
        CLz :: T = 0.4    # Dissipation parameter
        CLb :: T = 0.7    # Dissipation parameter
         Cτ :: T = 400.0  # Dissipation parameter

        CDe :: T = 2.0    # Dissipation parameter
        Ke₀ :: T = 1e-6   # Interior diffusivity for salinity

       CK_U :: T = 0.1    # Diffusivity parameter
       CK_T :: T = 0.1    # Diffusivity parameter
       CK_e :: T = 0.1    # Diffusivity parameter

        KU₀ :: T = 1e-6   # Interior viscosity for velocity
        KT₀ :: T = 1e-7   # Interior diffusivity for temperature
        KS₀ :: T = 1e-9   # Interior diffusivity for salinity
end

mutable struct State{T} <: FieldVector{6, T}
    Fu :: T
    Fv :: T
    Fθ :: T
    Fs :: T
    Fb :: T
     h :: T
    function State(T=Float64)
        new{T}(0, 0, 0, 0, 0, 0)
    end
end

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
    m.state.h  = ModularKPP.mixing_depth(m)
    return nothing
end

struct Model{KP, HP, SP, S, BC, TS, G, T} <: AbstractModularKPPModel{KP, HP, Nothing, TS, G, T}
      @add_clock_grid_timestepper
        solution :: S
             bcs :: BC
             tke :: TKEParameters{T}
     diffusivity :: KP
     mixingdepth :: HP
        kprofile :: SP
       constants :: Constants{T}
           state :: State{T}
end

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
             tke = TKEParameters(),
     diffusivity = ModularKPP.LMDDiffusivity(),
     mixingdepth = ModularKPP.LMDMixingDepth(),
        kprofile = ModularKPP.DiffusivityShape(),
         stepper = :ForwardEuler,
    )

    solution = Solution((CellField(grid) for i=1:nsol)...)

    bcs = (
        U = DefaultBoundaryConditions(eltype(grid)),
        V = DefaultBoundaryConditions(eltype(grid)),
        T = DefaultBoundaryConditions(eltype(grid)),
        S = DefaultBoundaryConditions(eltype(grid)),
        e = DefaultBoundaryConditions(eltype(grid))
    )

    Kϕ = (U=KU, V=KV, T=KT, S=KS, e=Ke)
    Rϕ = (U=RU, V=RV, T=RT, S=RS, e=Re)
    eq = Equation(K=Kϕ, R=Rϕ, update=update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, tke,
                    diffusivity, mixingdepth, kprofile, constants, State())
end

# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

@inline function mixing_time(m, i)
    return 400.0
    #N = sqrt( max(0, ∂B∂z(m, i)) )
    #@inbounds ũ = sqrt( max(0, m.solution.e[i]) ) + 1e-16
    #return @inbounds min(m.tke.Cτ, m.tke.CLz * abs(m.grid.zf[i]) / ũ, m.tke.CLb / N)
    #return @inbounds min(m.tke.Cτ, m.tke.CLb / N)
end

@inline K_mixing_time(m::AbstractModel{TS, G, T}, i) where {TS, G, T} =
    max(zero(T), mixing_time(m, i) * onface(m.solution.e, i)) + m.tke.KU₀

@inline oncell(f::Function, m, i) = 0.5 * (f(m, i) + f(m, i+1))
@inline onface(f::Function, m, i) = 0.5 * (f(m, i) + f(m, i-1))

@inline ∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α,
                            m.constants.β, i)

@inline production(m, i) = KU(m, i) * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2)

@inline buoyancy_flux(m, i) = - m.constants.g * (
      m.constants.α * KT(m, i) * ∂z(m.solution.T, i)
    - m.constants.β * KS(m, i) * ∂z(m.solution.S, i)
    )

@inline dissipation(m, i) = @inbounds m.tke.CDe * m.solution.e[i] / mixing_time(m, i)

#
# Equation entry
#

@inline K(m, i) = @inbounds mixing_time(m, i) * onface(m.solution.e, i)

#=
@inline KU(m, i) = m.tke.KU₀ #+ m.tke.CK_U * K(m, i)
@inline KT(m, i) = m.tke.KT₀ #+ m.tke.CK_T * K(m, i)
@inline KS(m, i) = m.tke.KS₀ #+ m.tke.CK_T * K(m, i)
@inline Ke(m, i) = m.tke.Ke₀ #+ m.tke.CK_e * K(m, i)
=#

@inline KU(m, i) = ModularKPP.KU(m, i)
@inline KT(m, i) = ModularKPP.KT(m, i)
@inline KS(m, i) = ModularKPP.KS(m, i)
@inline Ke(m, i) = ModularKPP.KU(m, i)

const KV = KU

@inline RU(m, i) = @inbounds   m.constants.f * m.solution.V[i]
@inline RV(m, i) = @inbounds - m.constants.f * m.solution.U[i]
@inline RT(m, i) = 0
@inline RS(m, i) = 0


@inline Re(m, i) = oncell(production, m, i) + oncell(buoyancy_flux, m, i) - dissipation(m, i)

end # module
