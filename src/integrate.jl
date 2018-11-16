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
  @. T += dt * insolation(z, shortwave, longwave, consts.λsw, consts.λw) # distribute solar insolation
  T[end] += dt * (latentflux+sensibleflux)/dz         # add latent and sensible flux
  S[end] += dt * S[end]*(evap-precip)/dz
  nothing
end

step_TS!(p, c, dt, forcingargs...) = step_TS!(p.T, p.S, p.z, p.dz, c, dt, forcingargs...)
step_TS!(m, dt, forcingargs...) = step_TS!(m.profile, m.constants, dt, forcingargs...)
 
insolation(z, shortwave, longwave, dsw, dlw) = shortwave*exp(z/dsw) + longwave*exp(z/dlw)

function step_uv!(u, v, z, dz, f, ρ₀, imix, dt, τˣ, τʸ)
  #=
  Some math:
  ∂t(e^ift*U) = e^ift * Gz/ρ₀
  => e^if(t+dt)U = dt * e^(ift) * Gz/ρ₀
  => U = dt * e^(-ifdt) * Gz/ρ₀
  => u = dt * (+cos(fdt)*stress_x + sin(fdt)*stress_y) / (dz*ρ₀)
  => v = dt * (-sin(fdt)*stress_x + cos(fdt)*stress_y) / (dz*ρ₀)
  =#

  # Forward euler, exponential integral step.
  h = -z[imix] # distribute momentum uniformly throughout mixed layer
  @views @. u[imix:end] += dt * (sin(f*dt)*τʸ + cos(f*dt)*τˣ) / (ρ₀*h)
  @views @. v[imix:end] += dt * (cos(f*dt)*τʸ + sin(f*dt)*τˣ) / (ρ₀*h)
  nothing
end

step_uv!(p, c, args...) = step_uv!(p.u, p.v, p.z, p.dz, c.f, c.ρ₀, args...)
step_uv!(model, args...) = step_uv!(model.profile, model.constants, args...)
