@inline oncell(f::Function, m, i) = (f(m, i) + f(m, i+1)) / 2
@inline onface(f::Function, m, i) = (f(m, i) + f(m, i-1)) / 2

@inline sqrt_e(m, i) = @inbounds maxsqrt(m.solution.e[i])

@inline ∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α,
                          m.constants.β, i)

Base.@kwdef struct SimpleMixingLength{T} <: AbstractParameters
    CLz :: T = 0.4    # Dissipation parameter
    CLb :: T = Inf    # Dissipation parameter
    CLΔ :: T = 0.1    # Dissipation parameter
end

@inline function mixing_length(m::Model{<:SimpleMixingLength}, i)
    Ls = - m.mixing_length.CLz * m.grid.zf[i]
    Ls = isnan(Ls) ? Inf : Ls

    LN = m.mixing_length.CLb * onface(sqrt_e, m, i) / maxsqrt(∂B∂z(m, i))
    LN = isnan(LN) ? Inf : LN

    L = min(Ls, LN)
    L = L == Inf ? 0.0 : L

    L = max(L, m.mixing_length.CLΔ * m.grid.Δf)

    return L
end
