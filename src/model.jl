mutable struct Model{T,AF,TFI}
    t::Float64
    step::Int
    imix::Int
    params::Parameters
    ocean::Ocean{T} 
    forcing::Forcing{AF}
    finterp::ForcingInterpolant{TFI}
end

function Model(; forcing=Forcing(), params=Parameters(), ocean=Ocean())
    Model(0.0, 0, ocean.nz, params, ocean, forcing, ForcingInterpolant(forcing))
end

function updatevars!(model)
    updatedensity!(model)
    @. model.ocean.u = real(model.ocean.U)
    @. model.ocean.v = imag(model.ocean.U)
    nothing
end

mixedlayerdepth(model) = mixedlayerdepth(model.ocean.zᶠ, model.imix)
updatedensity!(model) = updatedensity!(model.ocean, model.params)
