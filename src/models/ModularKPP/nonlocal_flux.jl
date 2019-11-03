#####
##### Counter gradient flux model proposed by Large et al (1994)
#####

Base.@kwdef struct LMDCounterGradientFlux{T} <: AbstractParameters
    CNL :: T = 6.33 # Mass flux proportionality constant
end

const CGModel = Model{K, <:LMDCounterGradientFlux} where K

update_nonlocal_flux!(model::CGModel) = nothing

instantiate_plume(::LMDCounterGradientFlux, grid) = (T=nothing, S=nothing, W²=nothing)

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
    Cbw :: T = 2.86
     Ce :: T = 0.4
    Cew :: T = 0.572
     Cα :: T = 1.0
    Cστ :: T = 2.2
    Cσb :: T = 1.32
end

instantiate_plume(::DiagnosticPlumeModel, grid) = 
    (T=CellField(grid), S=CellField(grid), W²=FaceField(grid))

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

@inline w_standard_dev(m, i) = @inbounds w_standard_dev(ωb(m), ωτ(m), m.nonlocalflux.Cστ, 
                                                        m.nonlocalflux.Cσb, -m.grid.zc[i] / m.state.h)

#####
##### Entrainment
#####

@inline entrainment(z, Δz, Ce, h) = Ce * (1 / (Δz - z) + 1 / (Δz + z + h))

@inline entrainment(z, model) = 
    @inbounds entrainment(z, Δc(model.grid, model.grid.N), model.nonlocalflux.Ce, model.state.h)
                                                      
#####
##### Plume boundary conditions
#####

function set_tracer_plume_bc!(ϕ̆, ϕ, Qϕ, Cα, model)
    n = ϕ̆.grid.N

    # Surface layer model: √w² Δϕ̆ = - C Qϕ, where Δϕ̆ is plume excess.
    @inbounds ϕ̆[n] = ϕ[n] - Cα * Qϕ / w_standard_dev(model, n)

    return nothing
end

function set_vertical_momentum_plume_bc!(W²)
    @inbounds W²[W².grid.N+1] = 0
    return nothing
end

#####
##### Plume equations
#####

#####
##### Buoyancy
#####

@inline plume_buoyancy_excess(i, grid, T̆, S̆, T, S, α, β, g) =
    @inbounds g * (α*(T̆[i] - T[i]) - β*(S̆[i] - S[i]))

@inline plume_buoyancy_excess(i, grid, model) = 
    plume_buoyancy_excess(i, grid, model.state.plume.T, model.state.plume.S, 
                                   model.solution.T, model.solution.S, 
                                   model.constants.α, model.constants.β, model.constants.g) 

function update_nonlocal_flux!(model::Model{K, <:DiagnosticPlumeModel}) where K

    set_tracer_plume_bc!(model.state.plume.T, model.solution.T,
                         model.state.Qθ, model.nonlocalflux.Cα, model)

    set_tracer_plume_bc!(model.state.plume.S, model.solution.S,
                         model.state.Qs, model.nonlocalflux.Cα, model)

    set_vertical_momentum_plume_bc!(model.state.plume.W²)

    integrate_plume_equations!(model.state.plume.T, model.state.plume.S, model.state.plume.W²,
                               model.solution.T, model.solution.S, model.grid, model)

    return nothing
end

function integrate_plume_equations!(T̆, S̆, W̆², T, S, grid, model)

    n = grid.N

    # Vertical momentum at the nᵗʰ cell interface, approximating excess buoyancy at
    # interface with excess at top cell center:
    ΔB̆ᵢ₊₁ = plume_buoyancy_excess(n, grid, model)
    @inbounds W̆²[n] = -model.nonlocalflux.Cbw * Δf(grid, n) * ΔB̆ᵢ₊₁

    # Integrate from surface cell `N-1` downwards
    for i in n-1 : -1 : 1
        if W̆²[i+1] <= 0 # plume is stopped
            @inbounds W̆²[i] = 0
            @inbounds T̆[i] = T[i]
            @inbounds S̆[i] = S[i]
        else # plume still lives

            # Integrate temperature and salinity
            @inbounds T̆[i] = T̆[i+1] + Δc(grid, i+1) * entrainment(grid.zc[i+1], model) * (T̆[i+1] - T[i+1])
            @inbounds S̆[i] = S̆[i+1] + Δc(grid, i+1) * entrainment(grid.zc[i+1], model) * (S̆[i+1] - S[i+1])

            # Integrate vertical momentum
            ΔB̆ᵢ₊₁ = onface(i+1, grid, plume_buoyancy_excess, model)
                                                             
            @inbounds W̆²[i] = W̆²[i+1] - Δf(grid, i+1) * ( 
                                  model.nonlocalflux.Cbw * ΔB̆ᵢ₊₁
                                - model.nonlocalflux.Cew * entrainment(grid.zf[i+1], model) * W̆²[i+1])
        end
    end

    return nothing
end

#####
##### Mass flux contribution to tracer budgets
#####

maxzero(ϕ::T) where T = max(zero(T), ϕ)

mass_flux(m::Model{K, <:DiagnosticPlumeModel}, i) where K = 
    @inbounds -m.nonlocalflux.Ca * sqrt(maxzero(m.state.plume.W²[i]))

@inline M_Φ(i, grid, Φ, model) = @inbounds mass_flux(model, i) * Φ[i]

# Use upwards-biased difference to effect upwind differencing for a downward-travelling plume:
∂z_explicit_nonlocal_flux_T(m::Model{K, <:DiagnosticPlumeModel}, i) where K =
    @inbounds ∂z⁺(i, m.grid, M_Φ, m.state.plume.T, m)

∂z_explicit_nonlocal_flux_S(m::Model{K, <:DiagnosticPlumeModel}, i) where K =
    @inbounds ∂z⁺(i, m.grid, M_Φ, m.state.plume.S, m)
