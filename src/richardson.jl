#=
mixed layer schematic:
bottom      mld   surface
  -H ~~~~~~ -h ~~~~~ 0
  [1]       [ih]   [end]
  stable profile: ρ[1] > ρ[2] > ...

  note dρdz[2] = ρ[2]-d
=#

unstable(profile) = any(profile.dρdz .> 0)

# Brute force method
function convect!(p::Profile)
  calcdρdz!(profile)
  dρdz = profile.dρdz

  while unstable(profile)

    mixed = false
    i = p.nz

    while !mixed
      if ρ[i-1] < ρ[i] # unstable
        ρ[i-1:i] = 0.5*(ρ[i-1]+ρ[i])
        calcdρdz!(profile)
        mixed = true
      else
        i -= 1
      end
    end

  end

  nothing
end




function bulkrichardson(profile, ih)
   h = profile.z[ih]
  Δρ = profile.ρ[ih] - profile.ρ[ih-1]

  Δu = profile.u[ih] - profile.u[ih-1]
  Δv = profile.v[ih] - profile.v[ih-1]
  ΔUsq = Δu^2 + Δv^2

  -g*Δρ*h / (ρ₀*ΔUsq) # should be positive
end

function gradrichardson!(profile)
  @views @. profile.dρdz[2:end] = @. @views profile.ρ[2:end] - profile.ρ[1:end-1] / profile.dz
  @views @. profile.dvdz[2:end] = @. @views profile.v[2:end] - profile.v[1:end-1] / profile.dz
  @views @. profile.dudz[2:end] = @. @views profile.u[2:end] - profile.u[1:end-1] / profile.dz

  @. profile.dUdz = sqrt(profile.dudz^2 + profile.dvdz^2)
  @. profile.Ri = g*profile.dρdz / (ρ₀*profile.dUdz^2)

  nothing
end

function calcdρdz!(profile)
  @views @. profile.dρdz[2:end] = @. @views profile.ρ[2:end] - profile.ρ[1:end-1]
  profile.dρdz[1] = profile.dρdz[2]
  nothing
end

# positive dρdz (inc z, inc ρ) => unstable
# dρdz = 

