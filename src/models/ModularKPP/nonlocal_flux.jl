Base.@kwdef struct LMDCounterGradientFlux{T} <: AbstractParameters
    CNL :: T = 6.33 # Mass flux proportionality constant
end

function ∂NLT∂z(m::Model{K, <:LMDCounterGradientFlux}, i) where K
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Fθ, d(m, i+1), d(m, i), Δf(m.grid, i), m)
end

function ∂NLS∂z(m::Model{K, <:LMDCounterGradientFlux}, i) where K
    KPP.∂NL∂z(m.nonlocalflux.CNL, m.state.Fs, d(m, i+1), d(m, i), Δf(m.grid, i), m)
end


Base.@kwdef struct BulkPlumeParameters{T} <: AbstractParameters
     Ce :: T = 0.4
     Cμ :: T = 0.15
     Cb :: T = 0.5
     Cm :: T = 0.3
     Cα :: T = 1.0
     Cσ :: T = 1.0
    Cσb :: T = 1.0
end


σw(ωb, ωτ, Cσ, Cσb, d) = Cσ * (ωτ^3 + Cσb * ωb^3 * d)^(1/3) * (1 - d)^(1/2)

entrainment(Ce, h, Δz, z) = Ce * (- 1 / (z + Δz) + 1 / (h + z + Δz))

function plume_buoyancy(plume_T, plume_S, T, S, α, β, g, i)
    @inbounds g*(α*(plume_T[i] - T[i]) - β*(plume_S[i] - S[i]))
end

