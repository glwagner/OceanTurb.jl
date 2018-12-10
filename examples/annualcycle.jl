#=
Reproduce Glenn's simple example of a seasonal mixed layer cycle forced by
constant winds.

TODO: neutralize correctly... we want to ensure that incoming + outgoing = 0.
=#
using Pkg; Pkg.activate("..")

using OceanBoundaryLayerModels
using Statistics: mean

const day = OceanBoundaryLayerModels.day
const year = OceanBoundaryLayerModels.year

# Functions to construct a seasonal cycle + daily cycle of heating and cooling
seasonalcycle(t, dec, june) = dec + (june - dec) * 0.5 * (sin(2π/year*t) + 1)
sunlight(t) = max(cos(2π*t/day + π), 0) # 0 at night, sinusoid in daylight
nighttime(t) = floor(1 - sunlight(t))   # 1 at night, 0 during day

function getshortwave(t, dec_amp, june_amp) 
  seasonalcycle(t, dec_amp, june_amp) * sunlight(t)
end

function getlongwave(t, dec_amp, june_amp)
  seasonalcycle(t, dec_amp, june_amp) * nighttime(t)
end

# Compute neutral cycle
function neutralize(fcn, tdata, dec_init, june_init)
  correction = 0
  @eval begin
    test_forcing = $fcn.(tdata, $dec_init, $june_init)
    correction = mean(test_forcing)
  end
  dec_init-correction, june_init-correction
end

ndata = 365 * 24 # 24 samples per day
tdata = collect(range(0, stop=year, length=ndata))

# June and December incoming shortwave and outgoing longwave max and min
# Fluxes are positive downward
 Fsw_dec0 = 700 
Fsw_june0 = 1400

 Flw_dec0 = -117
Flw_june0 = -144

Fsw_dec, Fsw_june = neutralize(:getshortwave, tdata, Fsw_dec0, Fsw_june0)
Flw_dec, Flw_june = neutralize(:getlongwave, tdata, Flw_dec0, Flw_june0)

println("initial shortwave (june, dec): $Fsw_june0, $Fsw_dec0")
println("initial longwave (june, dec): $Flw_june0, $Flw_dec0")

println("final shortwave (june, dec): $Fsw_june, $Fsw_dec")
println("final longwave (june, dec): $Flw_june, $Flw_dec")
