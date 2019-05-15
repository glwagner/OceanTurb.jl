using Pkg; Pkg.activate("..")

using PyPlotPlus

using 
  OceanBoundaryLayerModels,
  OceanBoundaryLayerModels.IdealizedForcing

hour = OceanBoundaryLayerModels.hour
day = OceanBoundaryLayerModels.day

t = range(0, step=0.5hour, stop=4day)

shortwave_amp = 800
longwave_amp = -200

sw = shortwave_daily.(t, shortwave_amp)
lw = longwave_daily.(t, longwave_amp)

fig, ax = subplots()
plot(t/hour, sw, label="shortwave radiation")
plot(t/hour, lw, label="longwave radiation")

xlabel(L"t \, (\mathrm{hours})")
ylabel("forcing (\$ \\mathrm{W \\, m^{-2}} \$)")

legend()

cornerspines()
tightshow()
