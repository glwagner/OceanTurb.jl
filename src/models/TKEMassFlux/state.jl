mutable struct State{T, D, L, H, P}
       Qu :: T
       Qv :: T
       Qθ :: T
       Qs :: T
       Qb :: T
        K :: D
        ℓ :: L
        h :: H # boundary layer depth
    plume :: P # plumes
end

function State(grid, mixing_length, boundary_layer_depth, nonlocal_flux; T=Float64)
    mixing_length = instantiate_mixing_length(mixing_length)
    boundary_layer_depth = 4.0 #instantiate_boundary_layer_depth(boundary_layer_depth)
    plume = instantiate_plume(nonlocal_flux, grid)
    K = CellField(grid)
    return State((zero(T) for i=1:5)..., K, mixing_length, boundary_layer_depth, plume)
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

    zero_out_negative_tke!(m.solution.e)
    update_near_wall_tke!(m)
    update_boundary_layer_depth!(m)

    update_mixing_length!(m)
    update_diffusivity!(m)
    update_nonlocal_flux!(m)

    return nothing
end

# Fallbacks
instantiate_plume(nonlocal_flux, grid) = nothing #(T=nothing, S=nothing, e=nothing, W²=nothing)
update_nonlocal_flux!(model) = nothing
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

function update_diffusivity!(m)
    for i in 1:m.grid.N
        @inbounds m.state.K[i] = diffusivity_mixing_length(m, i) * sqrt_e(m, i)
    end

    # Neumann condition on eddy diffusivity over boundary
    @inbounds m.state.K[0] = m.state.K[1]
    @inbounds m.state.K[m.grid.N+1] = m.state.K[m.grid.N]

    return nothing
end
