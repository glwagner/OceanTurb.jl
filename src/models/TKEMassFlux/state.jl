mutable struct State{T} <: FieldVector{6, T}
    Qu :: T
    Qv :: T
    Qθ :: T
    Qs :: T
    Qb :: T
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

