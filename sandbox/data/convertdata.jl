using 
  NetCDF,
  EarthConstants,
  JLD2

filename = "SO_met_100day.nc"

ncinfo(filename)

t = ncread(filename, "time") * day
t_hrs = t / hour

# Radiation
shortwave = ncread(filename, "sw")
longwave = ncread(filename, "lw")

latentheating = ncread(filename, "qlat")
sensibleheating = ncread(filename, "qsens")

stress_x = ncread(filename, "tx")
stress_y = ncread(filename, "ty")

precip = ncread(filename, "precip")

filename = "example_southern_ocean_forcing.jld2"
file = jldopen(filename, "a")

file["t"]               = t
file["shortwave"]       = shortwave
file["longwave"]        = longwave
file["stress_x"]        = stress_x
file["stress_y"]        = stress_y
file["latentheating"]   = latentheating
file["sensibleheating"] = sensibleheating
file["precip"]          = precip

close(file)

#=
# Make plots
zeroline(t) = plot(t, 0t, color="k", alpha=0.2, linewidth=1)

close("all")
fig, axs = subplots(nrows=3, sharex=true, figsize=(12, 12))

sca(axs[1])
plot(t_hrs, stress_east;  label="East")
plot(t_hrs, stress_north; label="West")
zeroline(t_hrs)

xlabel("Time (hours)", labelpad=12)
ylabel("Wind stress (N/m\$^2\$)", labelpad=12)

cornerspines(horizontal="top")
legend(fontsize=12, labelspacing=0.5)

sca(axs[2])
plot(t_hrs, heating_latent;   label="Latent heat flux")
plot(t_hrs, heating_sensible; label="Sensible heat flux")
plot(t_hrs, shortwave;        label="Short wave radiation")
plot(t_hrs, longwave;         label="Long wave radiation")
zeroline(t_hrs)

ylabel("Heat flux (W/m)", labelpad=24)
legend(fontsize=10, labelspacing=0.2)

sidespine()

sca(axs[3])
plot(t_hrs, 1e6*precip)

xlabel("Time (hours)")
ylabel("Precipitation (\$\\mu\$m/s)", labelpad=24)

cornerspines()

tightshow()
=#
