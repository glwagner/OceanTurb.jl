mutable struct State{T, H, P}
        Qu :: T
        Qv :: T
        Qθ :: T
        Qs :: T
        Qb :: T
         h :: T
    h_crit :: H
     plume :: P
end

h_criterion(args...) = nothing

function State(diffusivity, nonlocalflux, mixingdepth, grid, T=Float64)
    plume = instantiate_plume(nonlocalflux, grid)
    h_crit = h_criterion(mixingdepth, grid)
    return State(zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
                 h_crit, plume)
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

    update_mixing_depth!(m)

    update_nonlocal_flux!(m)

    return nothing
end
