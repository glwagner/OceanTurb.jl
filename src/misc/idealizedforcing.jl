module IdealizedForcing

export
  shortwave_seasonal,
  longwave_seasonal,
  shortwave_daily,
  longwave_daily

using OceanTurb
using Statistics: mean

const day = OceanTurb.day
const year = OceanTurb.year

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

end # module
