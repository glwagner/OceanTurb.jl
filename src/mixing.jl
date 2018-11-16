#=
bottom      mld   surface
  -H ~~~~~~ -h ~~~~~ 0
  [1]       [ih]   [end]
  stable profile: ρ[1] > ρ[2] > ...

  note dρdz[2] = ρ[2]-d
=#

function pressenter()
  println("\nPress enter to continue.")
  chomp(readline())
end

unstable(profile) = any(s -> s>0, profile.dρdz)

#=
function convect!(prof::Profile)
  dz!(prof, :ρ)
  dρdz = prof.dρdz
  ρ = prof.ρ

  imix = nothing
  while unstable(prof)
    imix = _convect!(prof)
    dz!(prof, :ρ)
  end

  imix # exits with dρdz calculated.
end
=#

function convect!(prof, imix=prof.nz)
  avgρ = prof.ρ[imix] # mixed layer density
  while imix > 1 && avgρ >= prof.ρ[imix-1]
    imix -= 1 # quest downward
    avgρ = (prof.ρ[imix] + (prof.nz-imix)*avgρ) / (prof.nz-imix+1) # running average
  end
  homogenize!(prof, imix) # mix the profile from i to nz.
  imix
end

function homogenize!(prof, ih)
  @views prof.u[ih:prof.nz] .= mean(prof.u[ih:prof.nz])
  @views prof.v[ih:prof.nz] .= mean(prof.v[ih:prof.nz])
  @views prof.T[ih:prof.nz] .= mean(prof.T[ih:prof.nz])
  @views prof.S[ih:prof.nz] .= mean(prof.S[ih:prof.nz])
  @views prof.ρ[ih:prof.nz] .= mean(prof.ρ[ih:prof.nz])
  nothing
end

function bulk_richardson(prof, ih)
   h = -prof.z[ih]
  Δρ = prof.ρ[ih] - prof.ρ[ih-1] 
  Δu = prof.u[ih] - prof.u[ih-1]
  Δv = prof.v[ih] - prof.v[ih-1]
  ΔUsq = Δu^2 + Δv^2
  -g*Δρ*h / (ρ₀*ΔUsq) # should be positive...
end

function grad_richardson!(prof)
  ρ!(prof)
  dz!(prof, :ρ)
  dz!(prof, :u)
  dz!(prof, :v)
  noflux!(prof)
  @. prof.dUdz = sqrt(prof.dudz^2 + prof.dvdz^2)
  @. prof.Ri = g*prof.dρdz / (ρ₀*prof.dUdz^2)
  nothing
end

function dz!(prof, ss::Symbol)
  nz = prof.nz

  s = getfield(prof, ss)

  sdsdz = Symbol(:d, ss, :dz)
  dsdz = getfield(prof, sdsdz)

  @views @. dsdz[2:end] = s[2:end] - s[1:end-1]
  dsdz[1] = 0
  nothing
end

# positive dρdz (inc z, inc ρ) => unstable

function diffuse!(u, κ)
  @. u *= (1-2κ)
  u[1] = u[2] # no flux bottom boundary condition
  @views @. u[2:end-1] = κ*(u[1:end-2] + u[2:end])
  nothing
end

function noflux!(p::Profile)
  p.dρdz[1] = 0
  p.dudz[1] = 0
  p.dvdz[1] = 0
  p.u[1] = p.u[2]
  p.v[1] = p.v[2]
  p.T[1] = p.T[2]
  p.S[1] = p.S[2]
  p.ρ[1] = p.ρ[2]
  nothing
end

bulkmixing!(profile, Ri) = nothing
