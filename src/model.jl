struct Profile{El,A}
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

struct Constants
  f::Float64
  dsw::Float64
  dlw::Float64
  Cᵖ::Float64
  ρ₀::Float64
end

function Constants(; # verbose names for API:            
                                 latitude = 45,   # degrees
                     shortwavepenetration = 0.6,  # meters
                      longwavepenetration = 20,   # meters
                             heatcapacity = 3900, # ?
                               refdensity = 1025, # kg/m^3
                         # short names
                         dsw=shortwavepenetration, dlw=longwavepenetration, 
                         Cᵖ=heatcapacity, ρ₀=refdensity)
  latradians = latitude * π/180
  f = 2Ω*sin(latradians)
  Constants(f, dsw, dlw, Cᵖ, ρ₀)
end

struct Model{T,AP,AF,TFI}
  constants::Constants
  profile::Profile{T,AP} 
  forcing::Forcing{AF}
  finterp::ForcingInterpolant{TFI}
end

function Model(; forcing=Forcing(), constants=Constants(), H=100, nz=100, T=Float64)
  Model(constants, Profile(T, H, nz), forcing, ForcingInterpolant(forcing))
end
