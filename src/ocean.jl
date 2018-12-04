"""
    Ocean(T, H, nz)

Construct a 1D `Ocean`.
"""
struct Ocean{T}
  nz::Int
  dz::T
  z::AbstractArray{T}
  zF::AbstractArray{T}
  U::AbstractArray{Complex{T}}     
  u::AbstractArray{T}
  v::AbstractArray{T}     
  T::AbstractArray{T}     
  S::AbstractArray{T}     
  ρ::AbstractArray{T}  
  Ri::AbstractArray{T}
end

function Ocean(; Tel=Float64, H=400, nz=100)
  dz = H/nz
  zF = collect(Tel, range(-H, step=dz, stop=0)) # cell faces
  z  = 0.5*(zF[1:end-1] + zF[2:end]) # cell centers
  @zeros Tel (nz,) u v T S ρ Ri
  @zeros Complex{Tel} (nz,) U
  Ocean(nz, dz, z, zF, U, u, v, T, S, ρ, Ri)
end

Ocean(H, nz) = Ocean(Float64, H, nz)

function dz!(ocean::Ocean, ss::Symbol)
  nz = ocean.nz
  s = getfield(ocean, ss)
  sdsdz = Symbol(:d, ss, :dz)
  dsdz = getfield(ocean, sdsdz)
  @views @. dsdz[2:end] = s[2:end] - s[1:end-1]
  dsdz[1] = 0
  nothing
end

mixedlayerdepth(zF, imix) = -zF[imix] #+ 0.5*dz

density(T, S, ρ0, T0, S0, βT, βS) = ρ0*(1 - βT*(T-T0) + βS*(S-S0))
density(T, S, params) = density(T, S, params.ρ0, params.T0, params.S0, 
                                params.βT, params.βS)

function updatedensity!(ocean::Ocean, params)
  @. ocean.ρ = density(ocean.T, ocean.S,
                       params.ρ0, params.T0, params.S0, params.βT, params.βS)
  nothing
end

function noflux!(ocean::Ocean)
  ocean.U[1] = ocean.U[2]
  ocean.ρ[1] = ocean.ρ[2]
  nothing
end
