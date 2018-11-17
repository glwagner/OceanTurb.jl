function stepforward!(model, dt)
  # 1. Step forward temperature and salinity fields
  # 2. Resolve static instability through convective adjustment
  # 3. Step forward velocity fields
  # 4. Resolve mixed layer deepening due to turublent entrainment (bulk Richardson number criterion)
  # 5. ?
  nothing
end

"""
    step_TS!(model, dt)

Step forward the temperature and salinity fields.
"""
function step_TS!(T, S, z, dz, consts, dt, latentflux, sensibleflux, shortwave, longwave, evap, precip)
  @. T += dt * insolation(z, shortwave, longwave, consts.λswave, consts.λlwave) # distribute solar insolation
  T[end] += dt * (latentflux+sensibleflux)/dz         # add latent and sensible flux
  S[end] += dt * S[end]*(evap-precip)/dz
  nothing
end

step_TS!(p, c, dt, forcingargs...) = step_TS!(p.T, p.S, p.z, p.dz, c, dt, forcingargs...)
step_TS!(m, dt, forcingargs...) = step_TS!(m.profile, m.constants, dt, forcingargs...)
 
insolation(z, shortwave, longwave, λswave, λlwave) = shortwave*exp(z/λswave) + longwave*exp(z/λlwave)

"""
    step_U!(model, imix, dt, τˣ, τʸ)

Step forward `u` and `v` by `dt`, forcing the layer defined by `imix` by `τˣ` and `τʸ`.
"""
function step_U!(U, z, dz, f, ρ₀, imix, dt, τˣ, τʸ)
    #=
    Some math. With U = u + iv:
    ∂ₜ(e^{ift}*U) = e^{ift} * Gz/ρ₀
    => e^{if(t+dt)}*U = e^{ift}*U + dt * e^{ift} * Gz/ρ₀
    => U = e^{-ifdt}*U + dt * e^{-ifdt} * Gz/ρ₀
    =#

    # Forward euler, exponential integral step.
    @. U *= exp(-im*f*dt)

    # Forcing by surface stress, distributed over mixed layer depth `h`:
    h = -z[imix] # distribute momentum uniformly throughout mixed layer
    @views @. U[imix:end] += dt * exp(-im*f*dt)*(τˣ + im*τʸ) / (ρ₀*h)
    nothing
end

step_U!(p::Profile, c::Constants, imix, dt, τˣ, τʸ) = step_U!(p.U, p.z, p.dz, c.f, c.ρ₀, imix, dt, τˣ, τʸ)
step_U!(m::Model, imix, dt, τˣ, τʸ) = step_U!(m.profile, m.constants, imix, dt, τˣ, τʸ)
