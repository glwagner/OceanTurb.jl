Base.@kwdef struct LMDDiffusivity{T} <: AbstractParameters
     CKSL :: T = 0.1   # Surface layer fraction
       Cτ :: T = 0.4   # Von Karman constant

    Cstab :: T = 2.0   # Stable buoyancy flux parameter for wind-driven turbulence
    Cunst :: T = 6.4   # Unstable buoyancy flux parameter for wind-driven turbulence

       Cn :: T = 1.0   # Exponent for effect of stable buoyancy forcing on wind mixing
    Cmτ_U :: T = 0.25  # Exponent for effect of unstable buoyancy forcing on wind mixing of U
    Cmτ_T :: T = 0.5   # Exponent for effect of unstable buoyancy forcing on wind mixing of T
    Cmb_U :: T = 1/3   # Exponent for the effect of wind on convective mixing of U
    Cmb_T :: T = 1/3   # Exponent for effect of wind on convective mixing of T

     Cd_U :: T = 0.5   # Wind mixing regime threshold for momentum
     Cd_T :: T = 2.5   # Wind mixing regime threshold for tracers

     Cb_U :: T = 0.599 # Buoyancy flux parameter for convective turbulence
     Cb_T :: T = 1.36  # Buoyancy flux parameter for convective turbulence
    Cτb_U :: T = (Cτ / Cb_U)^(1/Cmb_U) * (1 + Cunst*Cd_U)^(Cmτ_U/Cmb_U) - Cd_U  # Wind stress parameter for convective turbulence
    Cτb_T :: T = (Cτ / Cb_T)^(1/Cmb_T) * (1 + Cunst*Cd_T)^(Cmτ_T/Cmb_T) - Cd_T  # Wind stress parameter for convective turbulence

      KU₀ :: T = 1e-6 # Interior viscosity for velocity
      KT₀ :: T = 1e-7 # Interior diffusivity for temperature
      KS₀ :: T = 1e-9 # Interior diffusivity for salinity
end

## ** The K-Profile-Parameterization **
K_KPP(h, 𝒲, d::T, p) where T = 0<d<1 ? max(zero(T), h * 𝒲 * shape(d, p)) : -zero(T)

𝒲_Holtslag(Cτ, Cτb, ωτ, ωb, d) = Cτ * (ωτ^3 + Cτb * d * ωb^3)^(1/3)
𝒲_Holtslag(m, i) = 𝒲_Holtslag(m.diffusivity.Cτ, m.diffusivity.Cτb, KPP.ωτ(m), KPP.ωb(m), KPP.d(m, i))

𝒲_LMD_unstable_U(m, i) = KPP.𝒲_unstable(
    m.diffusivity.CKSL, m.diffusivity.Cd_U,
    m.diffusivity.Cτ, m.diffusivity.Cunst,
    m.diffusivity.Cb_U, m.diffusivity.Cτb_U,
    m.diffusivity.Cmτ_U, m.diffusivity.Cmb_U,
    ωτ(m), ωb(m), d(m, i)
    )

𝒲_LMD_unstable_T(m, i) = KPP.𝒲_unstable(
    m.diffusivity.CKSL, m.diffusivity.Cd_T,
    m.diffusivity.Cτ, m.diffusivity.Cunst,
    m.diffusivity.Cb_T, m.diffusivity.Cτb_T,
    m.diffusivity.Cmτ_T, m.diffusivity.Cmb_T,
    ωτ(m), ωb(m), d(m, i)
    )

𝒲_LMD_stable(m, i) = KPP.𝒲_stable(
    m.diffusivity.Cτ, m.diffusivity.Cstab, m.diffusivity.Cn,
    ωτ(m), ωb(m), d(m, i)
    )

"Return the vertical velocity scale for momentum at face point i"
function 𝒲_LMD_U(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return 𝒲_LMD_unstable_U(m, i)
    else
        return 𝒲_LMD_stable(m, i)
    end
end

"Return the vertical velocity scale for tracers at face point i."
function 𝒲_LMD_T(m, i)
    if !isforced(m)
        return 0
    elseif isunstable(m)
        return 𝒲_LMD_unstable_T(m, i)
    else
        return 𝒲_LMD_stable(m, i)
    end
end

const 𝒲_LMD_V = 𝒲_LMD_U
const 𝒲_LMD_S = 𝒲_LMD_T

Base.@kwdef struct HoltslagDiffusivity{T} <: AbstractParameters
     Cτ :: T = 0.4
    Cτb :: T = 15.6
    KU₀ :: T = 1e-6 # Interior viscosity for velocity
    KT₀ :: T = 1e-7 # Interior diffusivity for temperature
    KS₀ :: T = 1e-9 # Interior diffusivity for salinity
end

