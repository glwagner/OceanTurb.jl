#####
##### Counter gradient flux model proposed by Large et al (1994)
#####

Base.@kwdef struct LMDCounterGradientFlux{T} <: AbstractParameters
    CNL :: T = 6.33 # Mass flux proportionality constant
end

const CGModel = Model{K, <:LMDCounterGradientFlux} where K

update_nonlocal_flux!(model::CGModel) = nothing

initialize_plumes(::LMDCounterGradientFlux, grid) = (T=nothing, S=nothing, W²=nothing)

mass_flux(m::CGModel, i) = 0

∂z_explicit_nonlocal_flux_T(m::Model{K, <:LMDCounterGradientFlux}, i) where K =
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Qθ, d(m, i+1), d(m, i), Δf(m.grid, i), m)

∂z_explicit_nonlocal_flux_S(m::Model{K, <:LMDCounterGradientFlux}, i) where K =
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Qs, d(m, i+1), d(m, i), Δf(m.grid, i), m)

#####
##### Diagnostic plume model
#####

Base.@kwdef struct DiagnosticPlumeModel{T} <: AbstractParameters
     Ca :: T = 0.1
     Cw :: T = 2.86
     Ce :: T = 0.4
    Cew :: T = 0.2
     Cα :: T = 1.0
    Cστ :: T = 2.2
    Cσb :: T = 1.32
end

initialize_plumes(::DiagnosticPlumeModel, grid) = 
    (T=CellField(grid), S=CellField(grid), W²=CellField(grid))

#####
##### Empirical standard deviation
#####

@inline function w_standard_dev(ωb, ωτ, Cστ, Cσb, d::T) where T
    if 0 < d < 1
        return (Cστ * ωτ^3 + Cσb * ωb^3 * d)^(1/3) * (1 - d)^(1/2)
    else
        return zero(T)
    end
end

@inline w_standard_dev(m, i) = w_standard_dev(ωb(m), ωτ(m), m.nonlocalflux.Cστ, 
                                              m.nonlocalflux.Cσb, d(m, i))

#####
##### Entrainment
#####

@inline entrainment(Ce, h, Δz, z) = Ce * (- 1 / (z + Δz) + 1 / (h + z + Δz))

@inline entrainment(i, grid, model) = @inbounds entrainment(model.nonlocalflux.Ce, model.state.h, 
                                                            Δc(grid, 1), grid.zc[i])

#####
##### Plume boundary conditions
#####

function set_tracer_plume_bc!(ϕ̆, ϕ, Qϕ, Cα, model)
    n = ϕ̆.grid.N

    # Surface layer model: √w² Δϕ̆ = - C Qϕ, where Δϕ̆ is plume excess.
    @inbounds ϕ̆[n] = ϕ[n] - Cα * Qϕ / w_standard_dev(model, n)

    # Set halo point just in case to ensure zero mass flux.
    @inbounds ϕ̆[n+1] = ϕ̆[n]
    return nothing
end

function set_vertical_momentum_plume_bc!(W²)
    @inbounds W²[W².grid.N+1] = 0
    @inbounds W²[W².grid.N] = 0
    return nothing
end

#####
##### Plume equations
#####

#####
##### Buoyancy
#####

plume_buoyancy_excess(i, grid, T̆, S̆, T, S, α, β, g) =
    @inbounds g * (α*(T̆[i] - T[i]) - β*(S̆[i] - S[i]))

function update_nonlocal_flux!(model::Model{K, <:DiagnosticPlumeModel}) where K

    set_tracer_plume_bc!(model.state.plumes.T, model.solution.T,
                         model.state.Qθ, model.nonlocalflux.Cα, model)

    set_tracer_plume_bc!(model.state.plumes.S, model.solution.S,
                         model.state.Qs, model.nonlocalflux.Cα, model)

    set_vertical_momentum_plume_bc!(model.state.plumes.W²)

    integrate_plume_equations!(model.state.plumes.T, model.state.plumes.S, model.state.plumes.W²,
                               model.solution.T, model.solution.S, model.grid, model)

    return nothing
end

function integrate_plume_equations!(T̆, S̆, W̆², T, S, grid, model)
    # Integrate from surface cell `N` downwards
    for i in grid.N-1:-1:1
        if i > 2 && W̆²[i+1] == 0 # Plume vanishes at i+1 or above: set plume quantities to zero below.
            @inbounds W̆²[i] = 0
            @inbounds T̆[i] = 0
            @inbounds S̆[i] = 0
        else # We have a plume
            @inbounds T̆[i] = T̆[i+1] + Δc(grid, i+1) * entrainment(i+1, grid, model) * (T̆[i+1] - T[i+1])
            @inbounds S̆[i] = S̆[i+1] + Δc(grid, i+1) * entrainment(i+1, grid, model) * (S̆[i+1] - S[i+1])

            ΔB̆ᵢ₊₁ = plume_buoyancy_excess(i+1, grid, T̆, S̆, T, S, model.constants.α, model.constants.β, 
                                          model.constants.g)
            @inbounds W̆²[i] = (W̆²[i+1] + model.nonlocalflux.Cw * Δc(grid, i+1) * 
                (ΔB̆ᵢ₊₁ - model.nonlocalflux.Cew * entrainment(i+1, grid, model) * W̆²[i+1]))

            if W̆²[i] < 0 # Plume energy is negative:
                W̆²[i] = 0 # Stop the plume
            end
        end
    end

    return nothing
end

#####
##### Mass flux contribution to tracer budgets
#####

mass_flux(m::Model{K, <:DiagnosticPlumeModel}, i) where K = 
    @inbounds -m.nonlocalflux.Ca * sqrt(m.state.plumes.W²[i])

@inline M_Φ(i, grid, Φ, model) = @inbounds mass_flux(model, i) * Φ[i]

# Use upwards-biased difference to effect upwind differencing for downward-travelling plumes:
∂z_explicit_nonlocal_flux_T(m::Model{K, <:DiagnosticPlumeModel}, i) where K =
    @inbounds ∂z⁺(i, m.grid, M_Φ, m.state.plumes.T, m)

∂z_explicit_nonlocal_flux_S(m::Model{K, <:DiagnosticPlumeModel}, i) where K =
    @inbounds ∂z⁺(i, m.grid, M_Φ, m.state.plumes.S, m)
