struct PriceWellerPinkelParameters
  critical_bulk_Ri::Float64
  critical_grad_Ri::Float64
  mix_grad_Ri::Float64
end




struct Parameters <: AbstractParameters
  f::Float64
  Cp::Float64
  ρ0::Float64
  g::Float64
  bulkRiᶜ::Float64
  gradRiᶜ::Float64
  gradRiᵐⁱˣ::Float64
  βᵀ::Float64 
  βˢ::Float64 
  T₀::Float64 
  S₀::Float64 
  dʳᵉᵈ::Float64    
  dᵇˡᵘᵉ::Float64
  fracʳᵉᵈ::Float64    
  fracᵇˡᵘᵉ::Float64
end


function Parameters(; # Verbose names for API:            
                                 latitude = 45,        # Latitude [degrees]
                             heatcapacity = 3990,      # Heat capacity of water [?]
                         referencedensity = 1.027e3,   # Reference density [g/m^3]
                                        g = 9.807,     # Gravitational acceleration [g/s^2]
                           criticalbulkRi = 0.65,      # Critical value for bulk Richardson number instability
                       criticalgradientRi = 0.25,      # Critical value for gradient Richardson number instability
                            gradientRimix = 0.3,       # Proportionality constant for gradient Ri mixing
                                       βᵀ = 1.67e-4,   # First thermal expansion coefficient [1/K]
                                       βˢ = 0.78e-3,   # Haline contraction coefficient [1/ppt]
                                       T₀ = 283,       # Reference temperature [K]
                                       S₀ = 35,        # Reference salinity [g/kg]
                     shortwavepenetration = 0.6,       # Penetration depth for long wave radiation [m]
                      longwavepenetration = 20,        # Penetration depth for long wave radiation [m]
                        shortwavefraction = 0.62,      # Penetration depth for long wave radiation [m]
                         longwavefraction = 0.38,      # Penetration depth for long wave radiation [m]
                         # Short names
                                 f = 2Ω*sin(latitude*π/180), 
                                Cᵖ = heatcapacity, 
                                ρ₀ = referencedensity,
                           bulkRiᶜ = criticalbulkRi, 
                           gradRiᶜ = criticalgradientRi,
                         gradRiᵐⁱˣ = gradientRimix,
                              dʳᵉᵈ = shortwavepenetration,   
                             dᵇˡᵘᵉ = longwavepenetration,   
                           fracʳᵉᵈ = shortwavefraction,   
                          fracᵇˡᵘᵉ = longwavefraction
                     )
                         
  Parameters(f, Cᵖ, ρ₀, g, bulkRiᶜ, gradRiᶜ, gradRiᵐⁱˣ, βᵀ, βˢ, T₀, S₀,
             dʳᵉᵈ, dᵇˡᵘᵉ, fracʳᵉᵈ, fracᵇˡᵘᵉ)
end

insolationprofile(z, frac_red, frac_blue, depth_red, depth_blue) = (
  frac_red*exp(z/depth_red) + frac_blue*exp(z/depth_blue))
