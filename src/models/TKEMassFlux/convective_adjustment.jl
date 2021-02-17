Base.@kwdef struct FluxProportionalConvectiveAdjustment{T} <: AbstractParameters
    Cᴬ :: T = 10.0
end

const ModelWithConvectiveAdjustment = Model{L, K, <:FluxProportionalConvectiveAdjustment} where {L, K}

@inline Cᴬ(m::ModelWithConvectiveAdjustment) = m.convective_adjustment.Cᴬ

@inline KU(m::ModelWithConvectiveAdjustment, i) = KU₀(m) + K_adjustment(Cᴷu(m, i), m, i)
@inline KV(m::ModelWithConvectiveAdjustment, i) = KV₀(m) + K_adjustment(Cᴷv(m, i), m, i)
@inline KT(m::ModelWithConvectiveAdjustment, i) = KT₀(m) + K_adjustment(CᴷT(m, i), m, i)
@inline KS(m::ModelWithConvectiveAdjustment, i) = KS₀(m) + K_adjustment(CᴷS(m, i), m, i)
@inline Ke(m::ModelWithConvectiveAdjustment, i) = Ke₀(m) + K_adjustment(Cᴷe(m, i), m, i)

@inline K_adjustment(Cᴷϕ, m, i) = 
    @inbounds ifelse(∂B∂z(m, i) < 0 && m.state.Qb > 0,
                     Cᴬ(m) * onface(m.solution.e, i)^2 / m.state.Qb,
                     Cᴷϕ * onface(m.state.K, i))
