struct Parameters
    f::Float64
    λˢʷ::Float64 # decay scale for shortwave absorption
    λˡʷ::Float64 # decay scale for longwave absorption
    Cᵖ::Float64
    ρ₀::Float64
    g::Float64
    bulkRiᶜ::Float64
    gradRiᶜ::Float64
    gradRiᵐⁱˣ::Float64
    βᵀ::Float64 
    βˢ::Float64 
    T₀::Float64 
    S₀::Float64 
    αᵥ::Float64 
end

function Parameters(; # Verbose names for API:            
                                 latitude = 45,        # Latitude [degrees]
                     shortwavepenetration = 0.6,       # Penetration depth for long wave radiation [m]
                      longwavepenetration = 20,        # Penetration depth for long wave radiation [m]
                             heatcapacity = 3900,      # Heat capacity of water [?]
                         referencedensity = 1.027e3,   # Reference density [g/m^3]
                                        g = 9.807,     # Gravitational acceleration [g/s^2]
                           criticalbulkRi = 0.65,      # Critical value for bulk Richardson number instability
                       criticalgradientRi = 0.25,      # Critical value for gradient Richardson number instability
                            gradientRimix = 0.3,       # Proportionality constant for gradient Ri mixing
                                       βᵀ = 1.67e-4,   # First thermal expansion coefficient [1/K]
                                       βˢ = 0.78e-3,   # Haline contraction coefficient [1/ppt]
                                       T₀ = 283,       # Reference temperature [K]
                                       S₀ = 35,        # Reference salinity [g/kg]
                                       αᵥ = 2.07e-4,   # Volumetric coefficient of thermal expansion for water [K⁻¹]
                         # Short names
                               λˢʷ = shortwavepenetration, 
                               λˡʷ = longwavepenetration,
                                 f = 2Ω*sin(latitude*π/180), 
                                Cᵖ = heatcapacity, 
                                ρ₀ = referencedensity,
                           bulkRiᶜ = criticalbulkRi, 
                           gradRiᶜ = criticalgradientRi,
                         gradRiᵐⁱˣ = gradientRimix
                     )
                         
    Parameters(f, λˢʷ, λˡʷ, Cᵖ, ρ₀, g, bulkRiᶜ, gradRiᶜ, gradRiᵐⁱˣ, βᵀ, βˢ, T₀, S₀, αᵥ)
end
