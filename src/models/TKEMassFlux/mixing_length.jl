Base.@kwdef struct SimpleMixingLength{T} <: AbstractParameters
    CLz :: T = 0.4
    CLb :: T = 0.1
    CLΔ :: T = 1.0
end

@inline function mixing_length_face(m::Model{<:SimpleMixingLength}, i)
    Ls = - m.mixing_length.CLz * m.grid.zf[i]

    LN = m.mixing_length.CLb * onface(sqrt_e, m, i) / maxsqrt(∂B∂z(m, i))
    LN = isnan(LN) ? Inf : LN

    L = min(Ls, LN)
    L = L == Inf ? 0.0 : L

    # Limit mixing length by some factor of the local cell width
    L = max(L, m.mixing_length.CLΔ * Δc(m.grid, i))

    return L
end

@inline function mixing_length_cell(m::Model{<:SimpleMixingLength}, i)
    Ls = - m.mixing_length.CLz * m.grid.zc[i]

    LN = m.mixing_length.CLb * sqrt_e(m, i) / oncell(sqrt_∂B∂z, m, i)
    LN = isnan(LN) ? Inf : LN

    L = min(Ls, LN)
    L = L == Inf ? 0.0 : L

    # Limit mixing length by some factor of the local face separation
    L = max(L, m.mixing_length.CLΔ * Δf(m.grid, i))

    return L
end
