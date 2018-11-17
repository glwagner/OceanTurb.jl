struct Profile{El,T}
    nz::Int
    dz::El
    z::AbstractArray{T}
    U::AbstractArray{Complex{T}}     
    u::AbstractArray{T}
    v::AbstractArray{T}     
    T::AbstractArray{T}     
    S::AbstractArray{T}     
    ρ::AbstractArray{T}  
    dUdz::AbstractArray{Complex{T}}
    dudz::AbstractArray{T}
    dvdz::AbstractArray{T}
    dρdz::AbstractArray{T}  
    Ri::AbstractArray{T}
end

function Profile(El, H, nz)
    dz = H/nz
    z = collect(El, range(-H+dz/2, step=dz, stop=-dz/2))
    @zeros El (nz,) u v T S ρ dudz dvdz dρdz Ri
    @zeros Complex{El} (nz,) U dUdz 
    Profile(nz, dz, z, U, u, v, T, S, ρ, dUdz, dudz, dvdz, dρdz, Ri)
end

Profile(H, nz) = Profile(Float64, H, nz)

struct Constants
    f::Float64
    λswave::Float64 # decay scale for shortwave absorption
    λlwave::Float64 # decay scale for longwave absorption
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
                         λswave=shortwavepenetration, λlwave=longwavepenetration,
                         f=2Ω*sin(latitude*π/180), Cᵖ=heatcapacity, ρ₀=refdensity)
                         
    Constants(f, λswave, λlwave, Cᵖ, ρ₀)
end

struct Model{T,AP,AF,TFI}
    constants::Constants
    profile::Profile{T,AP} 
    forcing::Forcing{AF}
    finterp::ForcingInterpolant{TFI}
end

function Model(; latitude=45, forcing=Forcing(), constants=Constants(latitude=latitude), H=400, nz=100, T=Float64)
    Model(constants, Profile(T, H, nz), forcing, ForcingInterpolant(forcing))
end

function dz!(prof, ss::Symbol)
    nz = prof.nz
    s = getfield(prof, ss)
    sdsdz = Symbol(:d, ss, :dz)
    dsdz = getfield(prof, sdsdz)
    @views @. dsdz[2:end] = s[2:end] - s[1:end-1]
    dsdz[1] = 0
    nothing
end

function updatevars!(model)
    @. model.profile.u = real(model.profile.U)
    @. model.profile.v = imag(model.profile.U)
    for var in (:U, :u, :v, :ρ)
        dz!(model.profile, var)
    end
    nothing
end
