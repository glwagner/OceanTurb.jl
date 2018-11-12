struct Profile{T}
  nz::Int
  z::T
  u::T
  v::T     
  T::T     
  S::T     
  ρ::T  
end

function Profile(Tel, H, nz)
  z = collect(Tel, range(-H, length=nz, stop=0))
  @zeros Tel (nz,) u v T S ρ
  Profile(nz, z, u, v, T, S, ρ)
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
