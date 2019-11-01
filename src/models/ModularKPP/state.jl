mutable struct State{T, H, U, W}
          Fu :: T
          Fv :: T
          Fθ :: T
          Fs :: T
          Fb :: T
           h :: T
      h_crit :: H
     plume_T :: U
     plume_S :: U
    plume_W² :: W
end

plumes(args...) = nothing, nothing, nothing
h_criterion(args...) = nothing

function State(diffusivity, nonlocalflux, mixingdepth, grid, T=Float64)
    plume_T, plume_S, plume_W² = plumes(nonlocalflux, grid)
    h_crit = h_criterion(mixingdepth, grid)
    return State(zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
                 h_crit, plume_T, plume_S, plume_W²)
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
    update_mixing_depth!(m)
    update_nonlocal_flux!(m)
    return nothing
end


