mutable struct Model{T,AF,TFI}
    t::Float64
    step::Int
    imix::Int
    Iˢʷ::AbstractArray{T}
    params::Parameters
    ocean::Ocean{T} 
    forcing::Forcing{AF}
    finterp::ForcingInterpolant{TFI}
end

function Model(; forcing=Forcing(), params=Parameters(), ocean=Ocean())
    Iˢʷ = insolationprofile(ocean, params) 
    Model(0.0, 0, ocean.nz, Iˢʷ, params, ocean, forcing, ForcingInterpolant(forcing))
end

function updatevars!(model)
    updatedensity!(model)
    @. model.ocean.u = real(model.ocean.U)
    @. model.ocean.v = imag(model.ocean.U)
    nothing
end

mixedlayerdepth(model) = mixedlayerdepth(model.ocean.zᴳ, model.imix)
updatedensity!(model) = updatedensity!(model.ocean, model.params)

function insolationprofile(ocean, params)
    @views @. (
        insolationprofile(ocean.zᴳ[2:end],   params.fracʳᵉᵈ, params.fracᵇˡᵘᵉ, params.dʳᵉᵈ, params.dᵇˡᵘᵉ)
        - insolationprofile(ocean.zᴳ[1:end-1], params.fracʳᵉᵈ, params.fracᵇˡᵘᵉ, params.dʳᵉᵈ, params.dᵇˡᵘᵉ))
end
