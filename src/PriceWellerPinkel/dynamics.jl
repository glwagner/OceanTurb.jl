mutable struct PriceWellerPinkelModel{T,AF,TFI} <: Model
  t::Float64
  step::Int
  imix::Int
  Iˢʷ::AbstractArray{T}
  params::Parameters
  ocean::Ocean{T} 
  forcing::Forcing{AF}
  finterp::ForcingInterpolant{TFI}
end






function PriceWellerPinkelModel(; forcing=Forcing(), params=Parameters(), 
                                  ocean=Ocean())
  Iˢʷ = insolationprofile(ocean, params) 
  PriceWellerPinkelModel(0.0, 0, ocean.nz, Iˢʷ, params, ocean, forcing, 
                         ForcingInterpolant(forcing))
end

function updatevars!(model)
  updatedensity!(model)
  @. model.ocean.u = real(model.ocean.U)
  @. model.ocean.v = imag(model.ocean.U)
  nothing
end

mixedlayerdepth(model) = mixedlayerdepth(model.ocean.zF, model.imix)
updatedensity!(model) = updatedensity!(model.ocean, model.params)

function insolationprofile(ocean, params)
  @views @. (
      insolationprofile(ocean.zF[2:end], params.frac_red, params.frac_blue, 
                        params.depth_red, params.depth_blue)
      - insolationprofile(ocean.zF[1:end-1], params.frac_red, params.frac_blue, 
                          params.depth_red, params.depth_blue))
                          
end


#=
bottom      mld   surface
  -H ~~~~~~ -h ~~~~~ 0
  [1]       [ih]   [end]
  stable profile: ρ[1] > ρ[2] > ...

  note dρdz[2] = ρ[2]-d
=#

function convect!(ocean::Ocean, imix)
  avgρ = ocean.ρ[imix] # mixed layer density
  while imix > 1 && avgρ >= ocean.ρ[imix-1]
    imix -= 1 # quest downward
    # Keep a running average of density
    avgρ = (ocean.ρ[imix] + (ocean.nz-imix)*avgρ) / (ocean.nz-imix+1)
  end
  homogenize!(ocean, imix) # mix the oile from i to nz.
  imix
end

function convect!(model::Model) 
  model.imix = convect!(model.ocean, model.imix)
  nothing
end

function homogenize!(ocean, i)
  @views ocean.U[i:end] .= mean(ocean.U[i:end])
  @views ocean.T[i:end] .= mean(ocean.T[i:end])
  @views ocean.S[i:end] .= mean(ocean.S[i:end])
  @views ocean.ρ[i:end] .= mean(ocean.ρ[i:end])
  nothing
end

gradRi(j, ρ, U, dz, g, ρ0) = -g*dz * (ρ[j]-ρ[j-1]) / (ρ0 * abs2(U[j]-U[j-1]))

gradRi(j, ocean, p) = gradRi(j, ocean.ρ, ocean.U, ocean.dz, p.g, p.ρ0)

function gradRi(j, ocean::Ocean, p::PriceWellerPinkelParameters) 
  gradRi(j, ocean.ρ, ocean.U, ocean.dz, p.g, p.ρ0)
end

function gradRi!(ocean::Ocean, p::Parameters, jrange)
  for j in jrange
    ocean.Ri[j] = j > 1 ? gradRi(j, ocean, p) : Inf
  end
  NaN2Inf!(ocean.Ri)
  nothing
end

function gradRi!(ocean::Ocean, params::Parameters, imix::Int)
  ocean.Ri .= Inf # Ri not defined at bottom point and in mixed layer
  gradRi!(ocean, params, 2:imix)
  nothing
end

gradRi!(model::Model) = gradRi!(model.ocean, model.params, model.imix)

function NaN2Inf!(a) 
  @. a[isnan(a)] = Inf
  nothing
end


"""
    shearmix!(u, j, Ri, Riᵐⁱˣ)

Mix values of `u` by the ratio `Ri/Riᵐⁱˣ` between `j` and `j-1`  to emulate 
shear mixing. `Ri` should always be less than 0.25; therefore `Riᵐⁱˣ < 0.25` 
implies imperfect mixing. The default value for `Riᵐⁱˣ` is 0.3.
"""
function shearmix!(u, j, Ri, Riᵐⁱˣ)
  uᵐⁱˣ = 0.5 * (1 - Ri[j]/Riᵐⁱˣ) * (u[j-1] - u[j])
  u[j-1] -= uᵐⁱˣ
  u[j]   += uᵐⁱˣ
  nothing
end

function shearmix!(ocean::Ocean, j, Riᵐⁱˣ)
  shearmix!(ocean.U, j, ocean.Ri, Riᵐⁱˣ)
  shearmix!(ocean.T, j, ocean.Ri, Riᵐⁱˣ)
  shearmix!(ocean.S, j, ocean.Ri, Riᵐⁱˣ)
  shearmix!(ocean.ρ, j, ocean.Ri, Riᵐⁱˣ)
  nothing
end

function gradientmix!(ocean, params, imix)
  gradRi!(ocean, params, imix)
  Riᵐⁱⁿ, jᵐⁱⁿ = findmin(ocean.Ri)

  while Riᵐⁱⁿ <= params.gradRiᶜ
    shearmix!(ocean, jᵐⁱⁿ, params.gradRiᵐⁱˣ)

    # Recompute Ri, taking care not to leave domain
    jlower = jᵐⁱⁿ > 1           ? jᵐⁱⁿ-1 : 1
    jupper = jᵐⁱⁿ+2 <= ocean.nz ? jᵐⁱⁿ+2 : ocean.nz
    gradRi!(ocean, params, jlower:jupper)

    Riᵐⁱⁿ, jᵐⁱⁿ = findmin(ocean.Ri) # recompute minimium Ri
  end

  nothing
end

gradientmix!(model::Model) = gradientmix!(model.ocean, model.params, model.imix)

function bulkRi(ρ, U, zF, dz, imix, g, ρ0)
  if imix > 1
    h = mixedlayerdepth(zF, imix)
    Δρ = ρ[imix] - ρ[imix-1] # negative if stably stratified
    ΔU = U[imix] - U[imix-1]
    return -g*Δρ*h / (ρ0*abs2(ΔU)) # > 0
  else
    return Inf
  end
end

bulkRi(o::Ocean, p::Parameters, imix) = bulkRi(o.ρ, o.U, o.zF, o.dz, imix, p.g, 
                                               p.ρ0)
bulkRi(m::Model) = bulkRi(m.ocean, m.params, m.imix)

function deepen!(ocean, imix)
  if imix > 1
    nmld = ocean.nz - imix + 1
    for fldname in (:U, :T, :S, :ρ)
      fld = getfield(ocean, fldname)
      deepen!(fld, imix, nmld)
    end
    imix -= 1
  end
  imix
end

function deepen!(fld, imix, nmld)
  newmixedlayerfld = (nmld*fld[imix] + fld[imix-1]) / (nmld+1)
  @views @. fld[imix-1:end] = newmixedlayerfld
  nothing
end

function bulkmix!(ocean, params, imix)
  Riᵇ = bulkRi(ocean, params, imix)
  while isfinite(Riᵇ) && Riᵇ <= params.bulkRiᶜ
    imix = deepen!(ocean, imix)
    Riᵇ = bulkRi(ocean, params, imix)
  end
  imix
end

function bulkmix!(m::Model) 
  m.imix = bulkmix!(m.ocean, m.params, m.imix)
  nothing
end
