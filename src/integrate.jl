function stepforward!(model, dt)

  # 1. Step forward temperature and salinity fields and update density profile
  FI = model.finterp # interpolator for forcing fields

  step_TS!(model, dt, 
            FI.latentflux(model.t),
            FI.sensibleflux(model.t),
            FI.shortwave(model.t),
            FI.longwave(model.t),
            FI.evap(model.t),
            FI.precip(model.t))

  updatedensity!(model) # Update model density profile from T-S

  # 2. Resolve static instability through convective adjustment and recalculate density
  convect!(model)

  # 3. Step forward velocity fields given surface wind stress
  step_U!(model, dt, FI.τˣ(model.t), FI.τʸ(model.t))

  # 4. Resolve mixed layer deepening due to turbulent entrainment (bulk Richardson number criterion)
  bulkmix!(model)

  # 5. Resolve turbulent mixing due to shear instability
  gradientmix!(model)

  nothing
end

function stepforward!(model, dt, nsteps::Int)
    for i = 1:nsteps
        stepforward!(model, dt)
    end
    nothing
end

"""
    step_TS!(model, dt)

Step forward the temperature and salinity fields.
"""
function step_TS!(T, S, z, dz, params, dt, latentflux, sensibleflux, shortwave, longwave, evap, precip)
  #@. T += dt * insolation(z, shortwave, longwave, params.λˢʷ, params.λˡʷ) # distribute solar insolation
  T[end] += dt * (latentflux+sensibleflux)/dz         # add latent and sensible flux
  S[end] += dt * S[end]*(evap-precip)/dz
  nothing
end

(step_TS!(ocean::Ocean, params::Parameters, dt, forcingargs...) = 
     step_TS!(ocean.T, ocean.S, ocean.z, ocean.dz, params, dt, forcingargs...))

step_TS!(model::Model, dt, forcingargs...) = step_TS!(model.ocean, model.params, dt, forcingargs...)

"""
    step_U!(model, dt, τˣ, τʸ)

Step forward `U = u + im*v` by `dt`, forcing the layer defined by `model.imix` by `τˣ + im*τʸ`.
"""
function step_U!(U, zᶠ, f, ρ₀, imix, τˣ, τʸ, dt) 
    #=
    Some math. With U = u + iv:
    ∂ₜ(e^{ift}*U) = e^{ift} * Gz/ρ₀
    => e^{if(t+dt)}*U₊₁ = e^{ift}*U + dt * e^{ift} * Gz/ρ₀
    => U = e^{-ifdt}*U + dt * e^{-ifdt} * Gz/ρ₀
    =#

    # Forward euler, exponential integral step.
    @. U *= exp(-im*f*dt)

    # Forcing by surface stress, distributed over mixed layer depth `h`:
    # Recall that ∂U = -Gz
    h = mixedlayerdepth(zᶠ, imix) # distribute momentum uniformly throughout mixed layer
    @views @. U[imix:end] -= dt * exp(-im*f*dt)*(τˣ + im*τʸ) / (ρ₀*h)
    nothing
end

step_U!(o::Ocean, p::Parameters, imix, τˣ, τʸ, dt) = step_U!(o.U, o.zᶠ, p.f, p.ρ₀, imix, τˣ, τʸ, dt)
step_U!(m::Model, dt, τˣ, τʸ) = step_U!(m.ocean, m.params, m.imix, τˣ, τʸ, dt)
