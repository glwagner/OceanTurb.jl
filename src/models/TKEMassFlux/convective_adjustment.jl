abstract type AbstractConvectiveAdjustmentParameters <: AbstractParameters end

const ModelWithConvectiveAdjustment = Model{L, K, <:AbstractConvectiveAdjustmentParameters} where {L, K}

@inline Cᴬu(m) = Cᴬu(m.convective_adjustment)
@inline Cᴬc(m) = Cᴬc(m.convective_adjustment)
@inline Cᴬe(m) = Cᴬe(m.convective_adjustment)

Base.@kwdef struct UnityPrandtlConvectiveAdjustment{T} <: AbstractConvectiveAdjustmentParameters  
    Cᴬ :: T = 1.0
end

@inline Cᴬu(p::UnityPrandtlConvectiveAdjustment) = p.Cᴬ
@inline Cᴬc(p::UnityPrandtlConvectiveAdjustment) = p.Cᴬ
@inline Cᴬe(p::UnityPrandtlConvectiveAdjustment) = p.Cᴬ

Base.@kwdef struct VariablePrandtlConvectiveAdjustment{T} <: AbstractConvectiveAdjustmentParameters  
    Cᴬu :: T = 0.5
    Cᴬc :: T = 1.0
    Cᴬe :: T = 0.5
end

@inline Cᴬu(p::VariablePrandtlConvectiveAdjustment) = p.Cᴬu
@inline Cᴬc(p::VariablePrandtlConvectiveAdjustment) = p.Cᴬc
@inline Cᴬe(p::VariablePrandtlConvectiveAdjustment) = p.Cᴬe

@inline KU(m::ModelWithConvectiveAdjustment, i) = KU₀(m) + K_adjustment(Cᴷu(m, i), Cᴬu(m), m, i)
@inline KV(m::ModelWithConvectiveAdjustment, i) = KV₀(m) + K_adjustment(Cᴷv(m, i), Cᴬu(m), m, i)
@inline KT(m::ModelWithConvectiveAdjustment, i) = KT₀(m) + K_adjustment(CᴷT(m, i), Cᴬc(m), m, i)
@inline KS(m::ModelWithConvectiveAdjustment, i) = KS₀(m) + K_adjustment(CᴷS(m, i), Cᴬc(m), m, i)
@inline Ke(m::ModelWithConvectiveAdjustment, i) = Ke₀(m) + K_adjustment(Cᴷe(m, i), Cᴬe(m), m, i)

@inline K_adjustment(Cᴷϕ, Cᴬϕ, m, i) = 
    @inbounds ifelse(∂B∂z(m, i) < 0 && m.state.Qb > 0,
                     Cᴬϕ * onface(m.solution.e, i)^2 / m.state.Qb,
                     Cᴷϕ * onface(m.state.K, i))
