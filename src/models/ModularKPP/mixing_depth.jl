#####
##### CVMix mixing depth model
#####

Base.@kwdef struct LMDMixingDepth{T} <: AbstractParameters
     CSL :: T = 0.1   # Surface layer fraction
     CRi :: T = 0.3   # Critical bulk Richardson number
     CKE :: T = 4.32  # Unresolved turbulence parameter
    CKE₀ :: T = 1e-11 # Minimum unresolved turbulence kinetic energy
end

function update_mixing_depth!(m::Model{K, NL, <:LMDMixingDepth}) where {K, NL}
    m.state.h  = mixing_depth(m)
    return nothing
end

bulk_richardson_number(m::AbstractModel, i) = KPP.bulk_richardson_number(
    m.solution.U, m.solution.V, m.solution.T, m.solution.S,
    m.state.Qb, m.mixingdepth.CKE, m.mixingdepth.CKE₀, m.mixingdepth.CSL, m.constants.g,
    m.constants.α, m.constants.β, i)

"""
    mixing_depth(model)

Calculate the mixing depth 'h' for `model`.
"""
function mixing_depth(m)
    ih₁ = m.grid.N + 1 # start at top.
    Ri₁ = bulk_richardson_number(m, ih₁) # should be 0.

    # Descend through grid until Ri rises above critical value
    while ih₁ > 1 && Ri₁ < m.mixingdepth.CRi
        ih₁ -= 1 # descend
        Ri₁ = bulk_richardson_number(m, ih₁)
    end

    # Edge cases:
    # 1. Mixing depth is 0:
    if ih₁ == m.grid.N + 1
        z★ = m.grid.zf[ih₁]

    # 2. Mixing depth is whole domain because Ri is always less than CRi:
    elseif ih₁ == 1 && Ri₁ < m.mixingdepth.CRi
        z★ = m.grid.zf[ih₁]

    # 3. Ri is infinite somewhere inside the domain.
    elseif !isfinite(Ri₁)
        z★ = m.grid.zc[ih₁]

    # Main case: mixing depth is in the interior.
    else # Ri₁ > CRi
        ΔRi = bulk_richardson_number(m, ih₁+1) - Ri₁ # <0 linearly interpolate to find h.
        # x = x₀ + Δx * (y-y₀) / Δy
        z★ = m.grid.zf[ih₁] + Δf(m.grid, ih₁) * (m.mixingdepth.CRi - Ri₁) / ΔRi
    end

    -z★ < 0 && @warn "mixing depth $(-z★) is negative"

    return -z★ # "depth" is negative height.
end

#####
##### ROMS mixing depth model
#####

Base.@kwdef struct ROMSMixingDepth{T} <: AbstractParameters
     CSL :: T = 0.1  # Surface layer fraction
     CRi :: T = 0.3  # Critical bulk Richardson number
     CKE :: T = 5.07 # Minimum unresolved turbulence kinetic energy
     CEk :: T = 0.0  # Turbulent Ekman depth parameter
end

h_criterion(::ROMSMixingDepth, grid) = FaceField(grid)

function update_mixing_depth!(m::Model{K, NL, <:ROMSMixingDepth}) where {K, NL}
    mixing_depth_criterion!(m.state.h_crit, m)
    m.state.h = mixing_depth(m)
    return nothing
end

h_weight(h, CSL, zf, i) = @inbounds -zf[i] / (CSL*h - zf[i])
h_weight(m, i) = h_weight(m.state.h, m.mixingdepth.CSL, m.grid.zf, i)

function h_kernel(U, V, T, S, CRi, CEk, g, α, β, f, i)
    @inbounds ∂z(U, i)^2 + ∂z(V, i)^2 - ∂B∂z(T, S, g, α, β, i)/CRi - CEk*f^2
end

h_kernel(m, i) = h_kernel(m.solution.U, m.solution.V, m.solution.T, m.solution.S,
                            m.mixingdepth.CRi, m.mixingdepth.CEk,
                            m.constants.g, m.constants.α, m.constants.β, m.constants.f, i)

function unresolved_kinetic_energy(m, i)
    @inbounds unresolved_kinetic_energy(-m.grid.zf[i],
        ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α, m.constants.β, i),
        m.state.Qb, m.mixingdepth.CKE, 0, m.constants.g, m.constants.α, m.constants.β)
end

"Calculate the mixing depth criterion function by integrating from z=0 downwards."
function mixing_depth_criterion!(h_crit, m)
    @inbounds h_crit[m.grid.N+1] = 0

    for i = m.grid.N:-1:1
        @inbounds h_crit[i] = h_crit[i+1] + h_weight(m, i) * h_kernel(m, i) * Δc(m.grid, i)
    end

    for i in eachindex(h_crit)
        @inbounds h_crit[i] -= unresolved_kinetic_energy(m, i) / m.grid.zf[i]
    end

    return nothing
end

linear_interp(y★, x₀, y₀, Δx, Δy) = x₀ + Δx * (y★ - y₀) / Δy

function mixing_depth(m::Model{K, NL, <:ROMSMixingDepth}) where {K, NL}
    ih₁ = findprev(x -> x<=0, m.state.h_crit.data, m.grid.N)
    @inbounds begin
        if ih₁ === nothing # Mixing depth is entire grid
            z★ = m.grid.zf[1]
        elseif ih₁ == m.grid.N # Mixing depth at surface?
            z★ = ifelse(m.state.h_crit[ih₁]==0, m.grid.zf[m.grid.N], m.grid.zf[m.grid.N+1])
        else # linearly interpolate
            # x = x₀ + Δx * (y-y₀) / Δy
            z★ = linear_interp(0, m.grid.zf[ih₁], m.state.h_crit[ih₁], Δf(m.grid, ih₁),
                                m.state.h_crit[ih₁+1] - m.state.h_crit[ih₁])
        end
    end

    return -z★
end
