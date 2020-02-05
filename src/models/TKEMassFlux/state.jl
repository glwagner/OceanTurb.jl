mutable struct State{T, M, H}
    Qu :: T
    Qv :: T
    Qθ :: T
    Qs :: T
    Qb :: T
    mixing_length :: M
    h :: H # boundary layer depth
end

function State(mixing_length, boundary_layer_depth; T=Float64)
    mixing_length = instantiate_mixing_length(mixing_length)
    boundary_layer_depth = instantiate_boundary_layer_depth(boundary_layer_depth)
    return State((zero(T) for i=1:5)..., mixing_length, boundary_layer_depth)
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

    update_mixing_length!(m)
    update_boundary_layer_depth!(m)
    update_near_wall_tke!(m)

    zero_out_negative_tke!(m.solution.e)

    return nothing
end

# Fallbacks
update_near_wall_tke!(m) = nothing
instantiate_mixing_length(mixing_length_parameters) = nothing
instantiate_boundary_layer_depth(boundary_layer_depth_parameters) = nothing
update_mixing_length!(m) = nothing
update_boundary_layer_depth!(m) = nothing

@inline function zero_out_negative_tke!(e)
    for i in eachindex(e)
        @inbounds e[i] = ifelse(e[i] < 0, zero(eltype(e)), e[i])
    end
    return nothing
end
