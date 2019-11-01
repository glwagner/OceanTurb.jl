# Fallbacks
update_nonlocal_flux!(m) = nothing

initialize_plumes(args...) = (T=nothing, S=nothing, W²=nothing)

massflux(model, i) = 0

Base.@kwdef struct LMDCounterGradientFlux{T} <: AbstractParameters
    CNL :: T = 6.33 # Mass flux proportionality constant
end

∂z_explicit_nonlocal_flux_T(m::Model{K, <:LMDCounterGradientFlux}, i) where K =
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Fθ, d(m, i+1), d(m, i), Δf(m.grid, i), m)

∂z_explicit_nonlocal_flux_S(m::Model{K, <:LMDCounterGradientFlux}, i) where K =
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Fs, d(m, i+1), d(m, i), Δf(m.grid, i), m)

mass_flux(m::Model{K, <:LMDCounterGradientFlux}, i) where K = 0

Base.@kwdef struct BulkPlumeParameters{T} <: AbstractParameters
     Cw :: T = 2.86
     Ce :: T = 0.4
    Cew :: T = 0.2
     Cα :: T = 1.0
    Cστ :: T = 2.2
    Cσb :: T = 1.32
end

initialize_plumes(::BulkPlumeParameters, grid) = 
    (T=CellField(grid), S=CellField(grid), W²=CellField(grid))

#####
##### Empirical standard deviation
#####

w_standard_dev(ωb, ωτ, Cσ, Cσb, d) = (Cστ * ωτ^3 + Cσb * ωb^3 * d)^(1/3) * (1 - d)^(1/2)

w_standard_dev(m, i) = w_standard_dev(ωb(m), ωτ(m), m.nonlocalflux.Cσ, 
                                      m.nonlocalflux.Cσb, d(m, i))

#####
##### Entrainment
#####

entrainment(Ce, h, Δz, z) = Ce * (- 1 / (z + Δz) + 1 / (h + z + Δz))

entrainment(model, i) = @inbounds entrainment(model.nonlocalflux.Ce, model.state.h, 
                                              Δc(model.grid, 1), model.grid.zc[i])

#####
##### Buoyancy
#####

plume_buoyancy_excess(T̆, S̆, T, S, α, β, g, i) =
    @inbounds g * (α*(T̆[i] - T[i]) - β*(S̆[i] - S[i]))

#####
##### Plume boundary conditions
#####

function scalar_plume_boundary_condition!(ϕ̆, ϕ, Qϕ, Cα)
    n = model.grid.N

    # Surface layer model: √w² Δϕ̆ = - C Qϕ, where Δϕ̆ is plume excess.
    @inbounds ϕ̆[n] = ϕ[n] - Cα * Qϕ / w_standard_dev(m, n)
    return nothing
end

vertical_momentum_plume_boundary_condition!(Ŵ) = @inbounds Ŵ[model.grid.N] = 0

#####
##### Plume equations
#####

function update_nonlocal_flux!(model)

    scalar_plume_boundary_condition!(model.state.plumes.T, model.solution.T,
                                     model.state.QT, model.nonlocalflux.Cα)

    scalar_plume_boundary_condition!(model.state.plumes.S, model.solution.S,
                                     model.state.QS, model.nonlocalflux.Cα)

    vertical_momentum_plume_boundary_condition!(model.state.plumes.W²)

    T = model.solution.T
    S = model.solution.S

    T̆  = model.state.plumes.T
    S̆  = model.state.plumes.S
    W̆² = model.state.plumes.W²

    # Integrate from surface cell `N` downwards
    for i in model.grid.N-1:-1:1
        @inbounds T̆[i] = T̆[i+1] + Δc(grid, i+1) * entrainment(model, i+1) * (T̆[i+1] - T[i+1])
        @inbounds S̆[i] = S̆[i+1] + Δc(grid, i+1) * entrainment(model, i+1) * (S̆[i+1] - S[i+1])

        ΔB̆ᵢ₊₁ = plume_buoyancy_excess(T̆, S̆, T, S, model.constants.α, model.constants.β, 
                                      model.constants.g, i+1)

        @inbounds W̆²[i] = (W̆²[i+1] + model.nonlocalflux.Cw * Δc(grid, i+1) * 
            (ΔB̆ᵢ₊₁ - model.nonlocalflux.Cew * entrainment(model, i+1) * W̆²[i+1]))
    end

    # Fill halos to ensure zero mass flux in first cell?
    #@inbounds T̆[0] = T̆[1]
    #@inbounds S̆[0] = S̆[1]
    #@inbounds W̆²[0] = 0

    return nothing
end

mass_flux(m::Model{K, <:BulkPlumeParameters}, i) where K = 
    @inbounds m.nonlocalflux.Ca * sqrt(model.state.plumes.W̆²[i])

M_Φ(i, grid, Φ, model) = @inbounds mass_flux(model, i) * Φ[i]

# Use upwards-biased difference to effect upwind differencing for downward-travelling plumes:
∂z_explicit_nonlocal_flux_T(model::Model{K, <:BulkPlumeParameters}, i) where K =
    @inbounds ∂z⁺(i, grid, M_Φ, m.state.plumes.T, model)

∂z_explicit_nonlocal_flux_S(model::Model{K, <:BulkPlumeParameters}, i) where K =
    @inbounds ∂z⁺(i, grid, M_Φ, m.state.plumes.S, model)
