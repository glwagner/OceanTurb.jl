module IdealizedForcing

export
  shortwave_seasonal,
  longwave_seasonal,
  shortwave_daily,
  longwave_daily

using OceanBoundaryLayerModels
using Statistics: mean

const day = OceanBoundaryLayerModels.day
const year = OceanBoundaryLayerModels.year

#=
This module constructs idealized forcing for test cases.

We condider:

  * A seasonal cycle of solar heating and cooling
=#

# Functions to construct a seasonal cycle + daily cycle of heating and cooling
seasonalcycle(t, dec, june) = dec + (june - dec) * 0.5 * (sin(2π/year*t) + 1)
sunlight(t) = max(cos(2π*t/day + π), 0) # 0 at night, sinusoid in daylight
nighttime(t) = floor(1 - sunlight(t))   # 1 at night, 0 during day

shortwave_daily(t, amp) = amp*sunlight(t)
longwave_daily(t, amp) = amp*nighttime(t)

function shortwave_seasonal(t, dec_amp, june_amp) 
  seasonalcycle(t, dec_amp, june_amp) * sunlight(t)
end

function longwave_seasonal(t, dec_amp, june_amp)
  seasonalcycle(t, dec_amp, june_amp) * nighttime(t)
end


# Compute neutral shortwave heating and longwave cooling given initial guess
#=
function neutralize_forcing(tdata, dec_shortwave_guess, june_shortwave_guess,
                                   dec_longwave_guess, june_longwave_guess)
  test_shortwave = getshortwave.(tdata, dec_shortwave_guess, june_shortwave_guess)
  test_longwave = getlongwave.(tdata, dec_guess, june_guess)

  # compute correction

  dec_init-correction, june_init-correction
end
=#

#ndata = 365 * 24 # 24 samples per day
#tdata = collect(range(0, stop=year, length=ndata))

#=
# June and December incoming shortwave and outgoing longwave max and min
# Fluxes are positive downward
 Fsw_dec0 = 700 
Fsw_june0 = 1400

 Flw_dec0 = -117
Flw_june0 = -144

Fsw_dec, Fsw_june = neutralize(:getshortwave, tdata, Fsw_dec0, Fsw_june0)
Flw_dec, Flw_june = neutralize(:getlongwave, tdata, Flw_dec0, Flw_june0)
=#

end # module
