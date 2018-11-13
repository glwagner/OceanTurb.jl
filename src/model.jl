struct Profile{T}
  nz::Int
  dz::Float64
  z::T
  u::T
  v::T     
  T::T     
  S::T     
  ρ::T  
  dudz::T
  dvdz::T
  dUdz::T
  dρdz::T  
  Ri::T
end

function Profile(Tel, H, nz)
  z = collect(Tel, range(-H, length=nz, stop=0))
  dz = z[2]-z[1]
  @zeros Tel (nz,) u v T S ρ dudz dvdz dUdz dρdz Ri
  Profile(nz, Float64(dz), z, u, v, T, S, ρ, dudz, dvdz, dUdz, dρdz, Ri)
end

Profile(H, nz) = Profile(Float64, H, nz)

struct Model{T}
  profile::Profile{T} 
  forcing::Forcing{T}
end

function Model(forcing::Forcing{T}, H, nz=100) where T
  profile = Profile(eltype(T), H, nz)
  Model(profile, forcing)
end
