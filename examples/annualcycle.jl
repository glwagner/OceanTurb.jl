#=
Reproduce Glenn's simple example of a seasonal mixed layer cycle forced by
constant winds.
=#

using PriceWellerPinkel

year = PriceWellerPinkel.year

ndata = 365 * 24 # samples per day
tdata = collect(range(0, stop=year, length=ndata))

qˢʷ_summer = 1400
qˢʷ_winter = 700 
qˡʷ_summer = -144
qˡʷ_winter = -117

function forcingexample(t)
    nothing
end
