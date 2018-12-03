"""
    Ocean(T, H, nz)

Construct a 1D `Ocean`.
"""
struct Ocean{T}
  nz::Int
  dz::T
  z::AbstractArray{T}
  zᴳ::AbstractArray{T}
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
  zᴳ = collect(Tel, range(-H, step=dz, stop=0)) # cell faces
  z  = 0.5*(zᴳ[1:end-1] + zᴳ[2:end]) # cell centers
  @zeros Tel (nz,) u v T S ρ Ri
  @zeros Complex{Tel} (nz,) U
  Ocean(nz, dz, z, zᴳ, U, u, v, T, S, ρ, Ri)
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

mixedlayerdepth(zᴳ, imix) = -zᴳ[imix] #+ 0.5*dz

density(T, S, ρ₀=1.027e3, T₀=283, S₀=35, βᵀ=1.67e-4, βˢ=0.78e-3) = ρ₀*(1 - βᵀ*(T-T₀) + βˢ*(S-S₀))
density(T, S, params::AbstractParameters) = density(T, S, params.ρ₀, params.T₀, params.S₀, params.βᵀ, params.βˢ)

function updatedensity!(ocean::Ocean, params)
  @. ocean.ρ = density(ocean.T, ocean.S,
                       params.ρ₀, params.T₀, params.S₀, params.βᵀ, params.βˢ)
  nothing
end

function noflux!(ocean::Ocean)
  ocean.U[1] = ocean.U[2]
  ocean.ρ[1] = ocean.ρ[2]
  nothing
end
