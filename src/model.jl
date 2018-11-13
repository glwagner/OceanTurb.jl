struct Profile{El, A}
  nz::Int
  dz::El
  z::A
  u::A
  v::A     
  T::A     
  S::A     
  ρ::A  
  dudz::A
  dvdz::A
  dUdz::A
  dρdz::A  
  Ri::A
end

function Profile(El, H, nz)
  dz = H/nz
  z = collect(El, range(-H+dz/2, step=dz, stop=-dz/2))
  @zeros El (nz,) u v T S ρ dudz dvdz dUdz dρdz Ri
  Profile(nz, dz, z, u, v, T, S, ρ, dudz, dvdz, dUdz, dρdz, Ri)
end

Profile(H, nz) = Profile(Float64, H, nz)

struct Model{El, T}
  profile::Profile{El, T} 
  forcing::Forcing{T}
end

function Model(forcing::Forcing{T}, H, nz=100) where T
  profile = Profile(eltype(T), H, nz)
  Model(profile, forcing)
end
