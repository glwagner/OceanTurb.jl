struct BoundaryLayerConstants <: Constants
  f::Float64
  Cp::Float64
  ρ0::Float64
  g::Float64
  βT::Float64 
  βS::Float64 
  T0::Float64 
  S0::Float64 
  depth_red::Float64
  depth_blue::Float64
  frac_red::Float64    
  frac_blue::Float64
end

function BoundaryLayerConstants(; 
  # Verbose names for API:            
                latitude = 45,        # Latitude [degrees]
            heatcapacity = 3990,      # Heat capacity of water [?]
        referencedensity = 1.027e3,   # Reference density [g/m^3]
                       g = 9.807,     # Gravitational acceleration [g/s^2]
                      βT = 1.67e-4,   # First thermal expansion coefficient [1/K]
                      βS = 0.78e-3,   # Haline contraction coefficient [1/ppt]
                      T0 = 283,       # Reference temperature [K]
                      S0 = 35,        # Reference salinity [g/kg]
    penetrationdepth_red = 0.6,       # Penetration depth for 'red' fraction of shortwave radiation [m]
   penetrationdepth_blue = 20,        # Penetration depth for 'blue' fraction of shortwave radiation [m]
   shortwavefraction_red = 0.38,      # Fraction of shortwave radiation in 'red' band 
  shortwavefraction_blue = 0.62,      # Fraction of shortwave radiation in 'blue' band
  # Short names:
       f = 2Ω*sin(latitude*π/180), 
      Cp = heatcapacity, 
      ρ0 = referencedensity,
   depth_red = penetrationdepth_red,
  depth_blue = penetrationdepth_blue,
   frac_red = shortwavefraction_red,
  frac_blue = shortwavefraction_blue,
  )
                         
  BoundaryLayerConstants(f, Cp, ρ0, g, βT, βS, T0, S0, depth_red, depth_blue, 
                     frac_red, frac_blue)
end
