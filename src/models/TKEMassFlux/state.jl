mutable struct State{T} <: FieldVector{6, T}
    Qu :: T
    Qv :: T
    Qθ :: T
    Qs :: T
    Qb :: T
    function State(T=Float64)
        new{T}(0, 0, 0, 0, 0)
    end
end

"""
    update_state!(model)

Update the top flux conditions and mixing depth for `model`
and store in `model.state`.
"""
function update_state!(m)
    m.state.Qu = getbc(m, m.bcs.U.top)
    m.state.Qv = getbc(m, m.bcs.V.top)
    m.state.Qθ = getbc(m, m.bcs.T.top)
    m.state.Qs = getbc(m, m.bcs.S.top)
    m.state.Qb = m.constants.g * (m.constants.α * m.state.Qθ - m.constants.β * m.state.Qs)
    return nothing
end
